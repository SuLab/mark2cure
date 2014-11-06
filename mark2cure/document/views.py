from django.conf import settings
from django.template import RequestContext
from django.shortcuts import render_to_response, get_object_or_404, redirect
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.contrib.auth.models import User

from mark2cure.document.models import Document, View, Section
from mark2cure.common.models import Task

from mark2cure.document.forms import DocumentForm, AnnotationForm
from mark2cure.document.utils import generate_results, create_from_pubmed_id
from mark2cure.document.serializers import TopUserFromViewsSerializer, AnnotationSerializer

from rest_framework import viewsets, generics


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
                              { 'task': task,
                                'doc': doc,
                                'sections': sections,
                                'user_profile': user_profile,
                                'task_type': 'concept-recognition' },
                              context_instance=RequestContext(request))


@login_required
@require_http_methods(['POST'])
def identify_annotations_submit(request, task_id, doc_id, section_id):
    '''
      This is broken out because there can be many submissions per document
      We don't want to use these submission to direct the user to elsewhere in the app
    '''
    task = get_object_or_404(Task, pk=task_id)
    document = get_object_or_404(Document, pk=doc_id)
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

    for section in sections:
        setattr(section, 'words', section.resultwords(user))
        setattr(section, 'user_annotations', section.latest_annotations(user))

    '''
      1) It's a GM doc with GM annotations used to score
      2) It has community contributions (from this experiment) for context
      3) It's a novel document annotated by the worker
    '''
    #activity = Activity(user=user, document=doc, task_type='cr', experiment=settings.EXPERIMENT if user_profile.mturk else None)
    #previous_activities_available = Activity.objects.filter(document=doc, task_type='cr', experiment=settings.EXPERIMENT if user_profile.mturk else None).exclude(user=user, user__userprofile__ignore=True).exists()
    previous_activities_available = False

    # Can't use a Document as a Golden Master if no GM annotations exist
    # (TODO) resolve the gm stuff
    #if doc.has_golden() and user_profile.current_gm:
    if doc.has_golden():
        results = {}
        score, true_positives, false_positives, false_negatives = generate_results(doc, user)
        results['score'] = score
        results['true_positives'] = true_positives
        results['false_positives'] = false_positives
        results['false_negatives'] = false_negatives

        activity.submission_type = 'gm'
        activity.precsion = score[0]
        activity.recall = score[1]
        activity.f_score = score[2]
        activity.save()

        user_profile.current_gm = False
        user_profile.save()

        return render_to_response('document/concept-recognition-results-gold.jade',
            { 'task': task,
              'doc': doc,
              'user_profile': user_profile,
              'sections': sections,
              'results': results,
              'task_type': 'concept-recognition' },
            context_instance=RequestContext(request))


    elif previous_activities_available:
        activity.submission_type = 'cc'
        activity.save()

        #user_profile.current_gm = False
        #user_profile.save()

        return render_to_response('document/concept-recognition-results-community.jade',
            { 'task': task,
              'doc': doc,
              'user_profile': user_profile,
              'sections': sections,
              'task_type': 'concept-recognition' },
            context_instance=RequestContext(request))


    elif not previous_activities_available:
        #activity.submission_type = 'na'
        #activity.save()

        #user_profile.current_gm = False
        #user_profile.save()

        return render_to_response('document/concept-recognition-results-not-available.jade',
            { 'doc': doc,
              'user_profile': user_profile,
              'sections': sections,
              'task_type': 'concept-recognition' },
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


class TopUserViewSet(generics.ListAPIView):
    serializer_class = TopUserFromViewsSerializer

    def get_queryset(self):
        doc_id = self.kwargs['doc_id']
        top_users = View.objects.filter(task_type='cr', completed=True, section__document__id=doc_id, experiment=settings.EXPERIMENT if self.request.user.userprofile.mturk else None).exclude(user=self.request.user, user__userprofile__ignore=True).values('user')
        top_users = [dict(y) for y in set(tuple(x.items()) for x in top_users)]
        return top_users[:4]


class AnnotationViewSet(generics.ListAPIView):
    serializer_class = AnnotationSerializer

    def get_queryset(self):
        section_id = self.kwargs['section_id']
        user_id = self.kwargs['user_id']

        section = get_object_or_404(Section, pk=section_id)
        user = get_object_or_404(User, pk=user_id)
        annotations = section.latest_annotations(user=user)
        return annotations
