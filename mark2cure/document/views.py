from django.template import RequestContext
from django.shortcuts import render_to_response, get_object_or_404, redirect
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.db.models import Q
from django.contrib.auth.models import User
from django.contrib.auth import authenticate, login, logout
from django.conf import settings
from django.core.mail import send_mail
from django.utils import simplejson

from mark2cure.document.models import *
from mark2cure.document.forms import DocumentForm, AnnotationForm, RefuteForm, CommentForm
from mark2cure.document.utils import generate_results, create_from_pubmed_id
from mark2cure.common.utils import get_timezone_offset, get_mturk_account
from mark2cure.document.serializers import TopUserFromViewsSerializer, AnnotationSerializer, RelationshipTypeSerializer

from rest_framework import viewsets, generics

from copy import copy
import oauth2 as oauth
import json, itertools, logging
logger = logging.getLogger(__name__)


@login_required
def list(request, page_num=1):
    doc_list = Document.objects.all()
    paginator = Paginator(doc_list, 25)

    try:
        docs = paginator.page(page_num)
    except PageNotAnInteger:
        docs = paginator.page(1)
    except EmptyPage:
        docs = paginator.page(paginator.num_pages)

    return render_to_response('document/list.jade',
                              {'docs': docs},
                              context_instance=RequestContext(request))


'''
  Views for completing the Concept Recognition task
'''
@login_required
def identify_annotations(request, doc_id):
    # If they're attempting to view or work on the document
    doc = get_object_or_404(Document, pk=doc_id)
    sections = doc.available_sections()

    '''
      Technically we may want a user to do the same document multiple times,
      just means that during the community consensus we don't include their own reults
      to compare against
    '''

    # Make sure there is an abstract to annotate
    if sections.filter(kind="a").count() is 0:
        return redirect('mark2cure.document.views.validate_concepts', doc.pk)

    doc.create_views(request.user, 'cr')
    user_profile = request.user.userprofile

    return render_to_response('document/concept-recognition.jade',
                              { 'doc': doc,
                                'sections' : sections,
                                'user_profile' : user_profile,
                                'task_type': 'concept-recognition' },
                              context_instance=RequestContext(request))



@login_required
@require_http_methods(["POST"])
def identify_annotations_submit(request, doc_id, section_id):
    '''
      This is broken out because there can be many submissions per document
      We don't want to use these submission to direct the user to elsewhere in the app
    '''
    section = get_object_or_404(Section, pk=section_id)
    # Save this as not complete until they all complete
    view = section.update_view(request.user, 'cr', False)

    form = AnnotationForm(request.POST, view)
    if form.is_valid():
        ann = form.save(commit=False)

        ann.view = view
        ann.type = "disease"
        ann.user_agent = request.META['HTTP_USER_AGENT']
        ann.player_ip = request.META['REMOTE_ADDR']

        if request.user.profile.mturk:
          ann.experiment = settings.EXPERIMENT

        ann.save()
        return HttpResponse(200)
    return HttpResponse(500)


def identify_annotations_results(request, doc_id):
    '''
      After a document has been submitted, show the results and handle score keeping details
    '''
    doc = get_object_or_404(Document, pk=doc_id)
    sections = doc.available_sections()
    user_profile = request.user.userprofile

    if not doc.is_complete(request.user, user_profile, sections):
      return redirect('mark2cure.document.views.identify_annotations', doc.pk)

    for section in sections:
      setattr(section, "words", section.resultwords(request.user))

      if user_profile.mturk:
        setattr(section, "user_annotations", section.annotations(request.user.username, experiment = settings.EXPERIMENT))
      else:
        setattr(section, "user_annotations", section.annotations(request.user.username))

    '''
      1) It's a GM doc with GM annotations used to score
      2) It has community contributions (from this experiment) for context
      3) It's a novel document annotated by the worker
    '''
    if user_profile.mturk:
        activity, created = Activity.objects.get_or_create(user=request.user, document=doc, task_type = 'cr', experiment=settings.EXPERIMENT)

        activity.experiment = settings.EXPERIMENT
        previous_activities_available = Activity.objects.filter(document = doc, task_type='cr', experiment=settings.EXPERIMENT).exclude(user=request.user).exists()
    else:
        activity, created = Activity.objects.get_or_create(user=request.user, document=doc, task_type = 'cr', experiment=None)
        previous_activities_available = Activity.objects.filter(document = doc, task_type='cr').exclude(user=request.user).exists()


    if doc.is_golden():
        results = {}
        score, true_positives, false_positives, false_negatives = generate_results(doc, request.user)
        results['score'] = score
        results['true_positives'] = true_positives
        results['false_positives'] = false_positives
        results['false_negatives'] = false_negatives
        if score[2] == 1.0 or score[2] == 0.0:
          send_mail('[Mark2Cure #6] HIT completion',
                        '{0} scored {1} on document id {2}'.format(request.user.pk, score[2], doc.pk),
                         settings.SERVER_EMAIL,
                         [email[1] for email in settings.MANAGERS])


        activity.submission_type = 'gm'
        activity.precsion = score[0]
        activity.recall = score[1]
        activity.f_score = score[2]
        activity.save()

        return render_to_response('document/concept-recognition-results-gold.jade',
            { 'doc': doc,
              'user_profile' : user_profile,
              'sections' : sections,
              'results' : results,
              'task_type': 'concept-recognition' },
            context_instance=RequestContext(request))


    elif previous_activities_available:
        activity.submission_type = 'cc'
        activity.save()
        return render_to_response('document/concept-recognition-results-community.jade',
            { 'doc': doc,
              'user_profile' : user_profile,
              'sections' : sections,
              'task_type': 'concept-recognition' },
            context_instance=RequestContext(request))

    elif not previous_activities_available:
        activity.submission_type = 'na'
        activity.save()
        return render_to_response('document/concept-recognition-results-not-available.jade',
            { 'doc': doc,
              'user_profile' : user_profile,
              'sections' : sections,
              'task_type': 'concept-recognition' },
            context_instance=RequestContext(request))


@require_http_methods(["POST"])
@login_required
def refute_section(request,  doc_id, section_id):
    view = get_object_or_404(View, section__pk = section_id, user = request.user, completed = True, task_type="cr")

    form = RefuteForm(request.POST)
    if form.is_valid():
      refute = form.save(commit=False)
      refute.view = view
      refute.save()
      return HttpResponse("Success")

    return HttpResponse('Unauthorized', status=401)


@require_http_methods(["POST"])
@login_required
def comment_document(request,  doc_id):
    doc = get_object_or_404(Document, pk=doc_id)

    form = CommentForm(request.POST)
    if form.is_valid():
      refute = form.save(commit=False)
      refute.document = doc
      refute.user = request.user
      refute.save()

      send_mail('[Mark2Cure #6] Document comment',
                '{0} commented: {1} on document id {2}'.format(request.user.pk, refute.message, doc.pk),
                settings.SERVER_EMAIL,
                [email[1] for email in settings.MANAGERS])

      return HttpResponse("Success")

    return HttpResponse('Unauthorized', status=401)


'''
  Views for completing the Verify Concept task
'''
@require_http_methods(["POST"])
@login_required
def validate_concepts(request, doc_id):
    doc = get_object_or_404(Document, pk=doc_id)
    relationships = doc.get_conceptrelation_entries_to_validate()

    if len(relationships) is 0:
      return redirect('mark2cure.document.views.identify_annotations', doc.pk)

    return render_to_response('document/verify-relationships.jade',
        { 'doc': doc,
          'relationships' : relationships,
          'task_type': 'validate-concepts' },
        context_instance=RequestContext(request))


@login_required
@require_http_methods(["POST"])
def validate_concepts_submit(request, doc_id):
    doc = get_object_or_404(Document, pk=doc_id)
    validating_cr = get_object_or_404(ConceptRelationship, pk=request.POST.get('concept_relationship'))
    overview = get_object_or_404(Section, kind='o', document=doc)

    view, vc = View.objects.get_or_create(section = overview, user = request.user)
    ann, ac = Annotation.objects.get_or_create(kind = 'r', view = view)

    concept_relationship = ConceptRelationship(
        concept = validating_cr.concept,
        relationship = validating_cr.relationship,
        target = validating_cr.target,
        annotation = ann,
        validate = validating_cr,
        confidence = 0 if request.POST.get('vote') == 'false' else 1)
    concept_relationship.save()

    return HttpResponse(200)


'''
  Views for other task types [IN PROGRESS]
'''
@login_required
def identify_concepts(request, doc_id):
    # If they're attempting to view or work on the document
    doc = get_object_or_404(Document, pk=doc_id)
    concepts = doc.get_concepts_for_classification()
    return render_to_response('document/identify-concepts.jade',
        { 'doc': doc,
          'concepts': concepts[:3],
          'task_type': 'validate-concepts' },
        context_instance=RequestContext(request))


@login_required
@require_http_methods(["POST"])
def identify_concepts_submit(request, doc_id):
    '''
      This is broken out because there can be many submissions per document
      We don't want to use these submission to direct the user to elsewhere in the app
    '''

    doc = get_object_or_404(Document, pk=doc_id)
    overview = doc.section_set.filter(kind = 'o').first()

    subject_concept = get_object_or_404(Concept, concept_id=request.POST["c_one"])
    object_concept = get_object_or_404(Concept, concept_id=request.POST["c_two"])

    view, vc = View.objects.get_or_create(section = overview, user = request.user)
    ann, ac = Annotation.objects.get_or_create(view = view, kind = 'r')

    for r in request.POST["relation"].split(","):
      relationship_type = get_object_or_404(RelationshipType, pk=r)

      ConceptRelationship.objects.get_or_create(
          concept = subject_concept,
          relationship = relationship_type,
          target = object_concept,
          annotation = ann)


    return HttpResponse(200)


'''
  Utility views for general document controls
'''
@login_required
@require_http_methods(['POST'])
def submit(request, doc_id):
    '''
      If the user if submitting results for a document an document and sections
    '''
    doc = get_object_or_404(Document, pk=doc_id)
    task_type = request.POST.get('task_type')


    if task_type == 'concept-recognition':
        doc.update_views(request.user, 'cr', True)
        return redirect('mark2cure.document.views.identify_annotations_results', doc.pk)

    elif task_type == 'validate-concepts':
        return redirect('mark2cure.document.views.validate_concepts', doc.pk)
    elif task_type == 'identify-concepts':
        return redirect('mark2cure.document.views.identify_concepts', doc.pk)


    else:
        doc.update_views(request.user, 'cr', True)
        return redirect('mark2cure.document.views.identify_annotations_results', doc.pk)


@login_required
@require_http_methods(['POST'])
def next(request, doc_id):
    '''
      If the user if submitting results for a document an document and sections
    '''
    # doc = get_object_or_404(Document, pk=doc_id)
    task_type = request.POST.get('task_type')
    doc = Document.objects.get_random_document()

    if task_type == 'concept-recognition':
        return redirect('mark2cure.document.views.identify_annotations', doc.pk)

    elif task_type == 'validate-concepts':
        return redirect('mark2cure.document.views.validate_concepts', doc.pk)
    elif task_type == 'identify-concepts':
        return redirect('mark2cure.document.views.identify_concepts', doc.pk)
    else:
        return redirect('mark2cure.document.views.identify_annotations', doc.pk)


@login_required
def delete(request, doc_id):
    ###############
    # (TODO) MUST BE ADMIN
    ###############

    # doc = get_object_or_404(Document, pk=doc_id)
    # doc.delete()
    return redirect('/document/')


@login_required
@require_http_methods(["POST"])
def create(request):
    '''
      Takes the document_id from POST and directs the
      user to that pubmed document, downloading it if nessesary
    '''
    form = DocumentForm(request.POST)
    if form.is_valid():
      doc = create_from_pubmed_id( request.POST['document_id'] )
      return redirect('mark2cure.document.views.identify_annotations', doc.pk)


class TopUserViewSet(generics.ListAPIView):
    serializer_class = TopUserFromViewsSerializer

    def get_queryset(self):
        doc_id = self.kwargs['doc_id']
        section_id = self.kwargs['section_id']

        if self.request.user.userprofile.mturk:
            top_users = View.objects.filter(task_type = 'cr', completed = True, section__document__id = doc_id, experiment = settings.EXPERIMENT).exclude(user = self.request.user).values('user').distinct()
        else:
            top_users = View.objects.filter(task_type = 'cr', completed = True, section__document__id = doc_id).exclude(user = self.request.user).values('user').distinct()

        return top_users


class AnnotationViewSet(generics.ListAPIView):
    serializer_class = AnnotationSerializer

    def get_queryset(self):
        section_id = self.kwargs['section_id']
        user_id = self.kwargs['user_id']
        return Annotation.objects.filter(view__section__id = section_id, view__user__id = user_id).all()[:10]


class RelationshipTypeViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows users to be viewed or edited.
    """
    queryset = RelationshipType.objects.all()
    serializer_class = RelationshipTypeSerializer


