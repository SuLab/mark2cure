from django.conf import settings
from django.template import RequestContext
from django.shortcuts import render_to_response, get_object_or_404, redirect
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, HttpResponseServerError
from django.contrib.auth.models import User

from django.template.response import TemplateResponse

from mark2cure.common.models import Task
from mark2cure.common.formatter import bioc_writer, bioc_as_json, apply_bioc_annotations

from .models import Document, Section, Annotation
from .forms import AnnotationForm
from .utils import generate_results, select_best_opponent
from .serializers import AnnotationSerializer

from rest_framework import generics

from brabeion import badges
import os
import random


def read_pubtator_bioc(request, pubmed_id, format_type):
    # When fetching via pubmed, include no annotaitons
    doc = get_object_or_404(Document, document_id=pubmed_id)
    reader = doc.get_pubtator()

    writer = bioc_writer(request)
    writer.collection = reader.collection

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')


def read_pubmed_bioc(request, pubmed_id, format_type):
    # When fetching via pubmed, include no annotaitons
    writer = bioc_writer(request)
    doc = get_object_or_404(Document, document_id=pubmed_id)
    writer = doc.get_pubtator()

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')


'''
  Views for completing the Concept Recognition task
'''


@login_required
@require_http_methods(['POST'])
def identify_annotations_submit(request, task_pk, section_pk):
    '''
      This is broken out because there can be many submissions per document
      We don't want to use these submission to direct the user to elsewhere in the app
    '''
    task = get_object_or_404(Task, pk=task_pk)
    section = get_object_or_404(Section, pk=section_pk)

    user_quest_rel = task.userquestrelationship_set.filter(user=request.user, completed=False).first()
    view = user_quest_rel.views.filter(section=section, completed=False).first()

    if view:
        form = AnnotationForm(data=request.POST or None)

        if form.is_valid():
            ann = form.save(commit=False)
            ann.view = view
            ann.save()
            return HttpResponse(200)

    return HttpResponseServerError()


@login_required
def identify_annotations_results_bioc(request, task_pk, doc_pk, format_type):
    task = get_object_or_404(Task, pk=task_pk)
    document = task.documents.filter(pk=doc_pk).first()
    if not document:
        return HttpResponseServerError()

    '''
        Try to find an optimal opponete to pair the player
        against. If one isn't available or none meet the minimum
        requirements then just tell the player they've
        annotated a new document
    '''
    opponent = select_best_opponent(task, document, request.user)
    if not opponent:
        # No other work has ever been done on this apparently
        # so we reward the user and let them know they were
        # first via a different template / bonus points
        uqr = task.userquestrelationship_set.filter(user=request.user).first()

        # Did the new user provide at least 1 annotation?
        # (TODO) Did the new annotations differ from pubtator?
        # it would make more sense to do that as a valid check of
        # "contribution" effort
        if Annotation.objects.filter(view__userquestrelationship=uqr, view__section__document=document).count() > 0:
            request.user.profile.rating.add(score=1000, user=None, ip_address=os.urandom(7).encode('hex'))
            badges.possibly_award_badge('points_awarded', user=request.user)
            return HttpResponseServerError('points_awarded')

        return HttpResponseServerError('no_points_awarded')

    # BioC Writer Response that will serve all partner comparison information
    writer = document.as_writer()
    writer = apply_bioc_annotations(writer, opponent)

    # Other results exist if other people have at least viewed
    # the quest and we know other users have at least submitted
    # results for this particular document
    player_views = []
    opponent_views = []
    for section in document.available_sections():
        # If paired against a player who has completed the task multiple times
        # compare the to the first instance of the person completing that Document <==> Quest
        # while taking the latest version of the player's

        uqr = task.userquestrelationship_set.filter(user=request.user).first()
        player_view = uqr.views.filter(user=request.user, section=section, completed=True).first()

        quest_rel = task.userquestrelationship_set.filter(user=opponent).first()
        opponent_view = quest_rel.views.filter(section=section, completed=True).first()

        player_views.append(player_view)
        opponent_views.append(opponent_view)

        # Save who the player was paired against
        player_view.opponent = opponent_view
        player_view.save()

    results = generate_results(player_views, opponent_views)
    score = results[0][2] * 1000
    if score > 0:
        request.user.profile.rating.add(score=score, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge('points_awarded', user=request.user)

    writer.collection.put_infon('flatter', random.choice(settings.POSTIVE_FLATTER) if score > 500 else random.choice(settings.SUPPORT_FLATTER))
    writer.collection.put_infon('points', str(score))

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')


'''
  Utility views for general document controls
'''

@login_required
@require_http_methods(['POST'])
def submit(request, task_id, doc_id):
    '''
      If the user if submitting results for a document an document and sections
    '''
    task = get_object_or_404(Task, pk=task_id)
    doc = get_object_or_404(Document, pk=doc_id)

    task_type = request.POST.get('task_type')

    if task_type == 'concept-recognition':
        task.complete_views(doc, request.user)
        return redirect('document:results', task.pk, doc.pk)
    else:
        task.complete_views(doc, request.user)
        return redirect('document:results', task.pk, doc.pk)


class AnnotationViewSet(generics.ListAPIView):
    serializer_class = AnnotationSerializer

    def get_queryset(self):
        section_id = self.kwargs['section_id']
        user_id = self.kwargs['user_id']

        section = get_object_or_404(Section, pk=section_id)
        user = get_object_or_404(User, pk=user_id)
        annotations = section.latest_annotations(user=user)
        return annotations
