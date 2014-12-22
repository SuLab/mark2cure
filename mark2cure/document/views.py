from django.template import RequestContext
from django.shortcuts import render_to_response, get_object_or_404, redirect
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.contrib.auth.models import User

from mark2cure.document.models import Document, Section
from mark2cure.common.models import Task

from mark2cure.document.forms import AnnotationForm
from mark2cure.document.utils import generate_results
from mark2cure.document.serializers import AnnotationSerializer

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
    user = request.user
    user_profile = user.userprofile

    '''
      Technically we may want a user to do the same document multiple times,
      just means that during the community consensus we don't include their own reults
      to compare against
    '''
    user_profile.user_agent = request.META['HTTP_USER_AGENT']
    user_profile.player_ip = request.META['REMOTE_ADDR']
    user_profile.save()

    return render_to_response('document/concept-recognition.jade',
                              {'task': task,
                               'doc': doc,
                               'sections': sections,
                               'user_profile': user_profile,
                               'task_type': 'concept-recognition'},
                              context_instance=RequestContext(request))


@login_required
@require_http_methods(['POST'])
def identify_annotations_submit(request, task_id, doc_id, section_id):
    '''
      This is broken out because there can be many submissions per document
      We don't want to use these submission to direct the user to elsewhere in the app
    '''
    task = get_object_or_404(Task, pk=task_id)
    section = get_object_or_404(Section, pk=section_id)

    user_quest_rel_views = task.userquestrelationship_set.get(user=request.user).views
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

    user_quest_rel_views = task.userquestrelationship_set.get(user=user).views

    # (TODO) Validate the number of required views for this document, etc...
    if not user_quest_rel_views.filter(section__document=doc, completed=True).exists():
        return redirect('mark2cure.document.views.identify_annotations', task.pk, doc.pk)

    '''
      1) It's a GM doc with GM annotations used to score
      2) It has community contributions (from this experiment) for context,
         only 1 is presented
      3) It's a novel document annotated by the worker
    '''

    # Select all (**including uncompleted**) other user started quests
    # that have been completed
    others_quest_relationships = task.userquestrelationship_set.exclude(user=user)
    gm_user = User.objects.get(username='Doc_G-man')
    selected_user = gm_user

    # Gather users from completed documents
    # that may come from uncompleted quests
    previous_users = []
    query = others_quest_relationships.exclude(user=gm_user)
    for quest_relationship in query:
        if quest_relationship.views.filter(section__document=doc, completed=True).exists():
            previous_users.append(quest_relationship.user)

    # Other results exist if other people have at least viewed
    # the quest and we know other users have at least submitted
    # results for this particular document
    user_views = []
    gm_views = []
    ctx = {'task': task,
           'doc': doc,
           'user_profile': user_profile,
           'task_type': 'concept-recognition'}

    if others_quest_relationships.exists() and \
            others_quest_relationships.filter(user=gm_user).exists() and \
            others_quest_relationships.get(user=gm_user).views.filter(section__document=doc, completed=True).exists():

        # There is an "expert's" annotations (GM) so
        # show those as the partner's
        for section in sections:
            user_view = user_quest_rel_views.get(section=section, completed=True)
            gm_quest_rel = others_quest_relationships.get(user=gm_user)

            gm_view = gm_quest_rel.views.get(section=section, completed=True)

            user_views.append(user_view)
            gm_views.append(gm_view)
            setattr(section, 'words', section.resultwords(user_view, gm_view))

        ctx['sections'] = sections
        ctx['partner'] = selected_user
        return show_comparison_results(request, user_views, gm_views, ctx)

    elif others_quest_relationships.exists() and len(previous_users):

        # No expert around so select a previous user at random
        random.shuffle(previous_users)
        selected_user = previous_users[0]

        for section in sections:
            user_view = user_quest_rel_views.get(section=section, completed=True)
            user_quest_rel = others_quest_relationships.filter(user=selected_user).first()
            gm_view = user_quest_rel.views.get(section=section, completed=True)

            user_views.append(user_view)
            gm_views.append(gm_view)
            setattr(section, 'words', section.resultwords(user_view, gm_view))

        ctx['sections'] = sections
        ctx['partner'] = selected_user
        return show_comparison_results(request, user_views, gm_views, ctx)

    else:
        # No other work has ever been done on this apparently
        # so we reward the user and let them know they were
        # first via a different template / bonus points
        for section in sections:
            user_view = user_quest_rel_views.get(section=section, completed=True)
            setattr(section, 'words', section.resultwords(user_view, False))

        request.user.profile.rating.add(score=1000, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge('points_awarded', user=request.user)

        ctx['sections'] = sections
        return render_to_response('document/concept-recognition-results-not-available.jade',
               ctx,
               context_instance=RequestContext(request))


def show_comparison_results(request, user_views, gm_views, ctx):
    # Take views from whoever the partner was
    # and use those to calculate the score (and assign
    # / reward as appropriate
    results = generate_results(user_views, gm_views)
    score = results[0][2] * 1000
    if score > 0:
        request.user.profile.rating.add(score=score, user=None, ip_address=os.urandom(7).encode('hex'))
        badges.possibly_award_badge('points_awarded', user=request.user)

    ctx['results'] = results
    return render_to_response('document/concept-recognition-results-partner.jade',
           ctx,
           context_instance=RequestContext(request))


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
        return redirect('mark2cure.document.views.identify_annotations_results', task.pk, doc.pk)
    else:
        task.complete_views(doc, request.user)
        return redirect('mark2cure.document.views.identify_annotations_results', task.pk, doc.pk)


class AnnotationViewSet(generics.ListAPIView):
    serializer_class = AnnotationSerializer

    def get_queryset(self):
        section_id = self.kwargs['section_id']
        user_id = self.kwargs['user_id']

        section = get_object_or_404(Section, pk=section_id)
        user = get_object_or_404(User, pk=user_id)
        annotations = section.latest_annotations(user=user)
        return annotations
