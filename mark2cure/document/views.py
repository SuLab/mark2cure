from django.conf import settings
from django.template import RequestContext
from django.shortcuts import render_to_response, get_object_or_404, redirect
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.contrib.auth.models import User

from django.template.response import TemplateResponse

from mark2cure.common.models import Task


from .models import Document, Section, Annotation
from .forms import AnnotationForm
from .utils import generate_results, select_best_opponent
from .serializers import AnnotationSerializer

from rest_framework import generics

from brabeion import badges
import os
import random


'''
  Views for completing the Concept Recognition task
'''


@login_required
def identify_annotations(request, task_id, doc_id, treat_as_gm=False):
    # If they're attempting to view or work on the document
    task = get_object_or_404(Task, pk=task_id)
    doc = get_object_or_404(Document, pk=doc_id)

    sections = doc.available_sections()

    ctx = { 'task': task,
            'doc': doc,
            'sections': sections,
            'user_profile': request.user.profile,
            'task_type': 'concept-recognition'}
    return TemplateResponse(request, 'document/concept-recognition.jade', ctx)


@login_required
@require_http_methods(['POST'])
def identify_annotations_submit(request, task_id, doc_id, section_id):
    '''
      This is broken out because there can be many submissions per document
      We don't want to use these submission to direct the user to elsewhere in the app
    '''
    task = get_object_or_404(Task, pk=task_id)
    section = get_object_or_404(Section, pk=section_id)

    user_quest_rel = task.userquestrelationship_set.filter(user=request.user, completed=False).first()
    user_quest_rel_views = user_quest_rel.views
    view = user_quest_rel_views.filter(section=section, completed=False).first()

    if view:
        form = AnnotationForm(request.POST)
        if form.is_valid():
            ann = form.save(commit=False)
            ann.view = view
            ann.save()
            return HttpResponse(200)

    return HttpResponse(500)


@login_required
def identify_annotations_results(request, task_id, doc_id):
    '''
      After a document has been submitted, show the results and handle score keeping details
    '''
    task = get_object_or_404(Task, pk=task_id)
    doc = get_object_or_404(Document, pk=doc_id)

    sections = doc.available_sections()
    user = request.user
    user_profile = user.userprofile

    # only take the one that has views
    user_quest_rel = task.userquestrelationship_set.filter(user=user, completed=False).latest()
    user_quest_rel_views = user_quest_rel.views

    # (TODO) Validate the number of required views for this document, etc...
    if not user_quest_rel_views.filter(section__document=doc, completed=True).exists():
        return redirect('document:read', task.pk, doc.pk)

    # Other results exist if other people have at least viewed
    # the quest and we know other users have at least submitted
    # results for this particular document
    player_views = []
    opponent_views = []
    ctx = {'task': task,
           'doc': doc,
           'user_profile': user_profile,
           'task_type': 'concept-recognition'}

    '''
        Try to find an optimal opponete to pair the player
        against. If one isn't available or none meet the minimum
        requirements then just tell the player they've
        annotated a new document
    '''
    opponent = select_best_opponent(task, doc, user)
    if opponent:

        for section in sections:
            # If paired against a player who has completed the task multiple times
            # compare the to the first instance of the person completing that Document <==> Quest
            # while taking the latest version of the player's

            player_view = user_quest_rel_views.filter(section=section, completed=True).first()

            # if this is 0 try the other
            # originates from Views not having timestamps for created / updated
            #if Annotation.objects.filter(view=player_view).count() == 0:

            quest_rel = task.userquestrelationship_set.filter(user=opponent).first()
            opponent_view = quest_rel.views.filter(section=section, completed=True).first()

            player_views.append(player_view)
            opponent_views.append(opponent_view)
            setattr(section, 'words', section.resultwords(player_view, opponent_view))

            # Save who the player was paired against
            player_view.opponent = opponent_view
            player_view.save()

        ctx['sections'] = sections
        ctx['partner'] = opponent
        return show_comparison_results(request, player_views, opponent_views, ctx)

    else:
        # No other work has ever been done on this apparently
        # so we reward the user and let them know they were
        # first via a different template / bonus points
        total_anns = 0
        for section in sections:
            user_view = user_quest_rel_views.filter(section=section, completed=True).first()
            setattr(section, 'words', section.resultwords(user_view, False))
            total_anns += Annotation.objects.filter(view=user_view).count()

        contributed = total_anns > 0
        if contributed:
            request.user.profile.rating.add(score=1000, user=None, ip_address=os.urandom(7).encode('hex'))
            badges.possibly_award_badge('points_awarded', user=request.user)

        ctx['sections'] = sections
        ctx['contributed'] = contributed
        return TemplateResponse(request,
                'document/concept-recognition-results-not-available.jade',
                ctx)


def show_comparison_results(request, user_views, gm_views, ctx, log_score=False):
    # Take views from whoever the partner was
    # and use those to calculate the score (and assign
    # / reward as appropriate
    results = generate_results(user_views, gm_views)
    score = results[0][2] * 1000
    if score > 0:
        request.user.profile.rating.add(score=score, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge('points_awarded', user=request.user)

    ctx['flatter'] = random.choice(settings.POSTIVE_FLATTER) if score > 500 else random.choice(settings.SUPPORT_FLATTER)
    ctx['results'] = results
    return TemplateResponse(request,
            'document/concept-recognition-results-partner.jade',
            ctx)


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
