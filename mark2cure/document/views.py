'''
doc views (controllers)
'''

from django.template import RequestContext
from django.shortcuts import render_to_response
from django.shortcuts import get_object_or_404, redirect
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.db.models import Q
from django.contrib.auth.models import User
from django.contrib.auth import authenticate, login, logout

from mark2cure.document.models import *
from mark2cure.document.forms import DocumentForm, AnnotationForm
from mark2cure.document.utils import create_from_pubmed_id, check_validation_status
from mark2cure.common.utils import get_timezone_offset, get_mturk_account


from rest_framework import viewsets
from mark2cure.document.serializers import RelationshipTypeSerializer

from copy import copy

import oauth2 as oauth
import json, itertools



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
def identify_annotations(request, doc_id):
    # If they're attempting to view or work on the document
    doc = get_object_or_404(Document, pk=doc_id)

    assignment_id = request.GET.get('assignmentId') #ASSIGNMENT_ID_NOT_AVAILABLE
    worker_id = request.GET.get('workerId')

    if doc.is_complete(request.user):
      return redirect('/document/'+ str(doc.pk) + '/results/' )

    if len(doc.section_set.filter(kind="a")) is 0:
      return redirect('/document/'+ str(doc.pk) + '/concepts/validate/' )

    # If mTurk user not logged in, make a new account for them and set the session
    if assignment_id == 'ASSIGNMENT_ID_NOT_AVAILABLE':
      logout(request)

    if worker_id and not request.user.is_authenticated():
      # If it's accepted and a worker that doesn't have an account
      user = get_mturk_account(worker_id)
      user = authenticate(username=user.username, password='')
      login(request, user)

    if request.user.is_authenticated():
      doc.update_views(request.user, 'cr')

    return render_to_response('document/concept-recognition.jade',
                              { 'doc': doc,
                                'task_type': 'concept-recognition',
                                'instruct_bool': 'block' if assignment_id == 'ASSIGNMENT_ID_NOT_AVAILABLE' else 'none',
                                'assignmentId': assignment_id},
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
          ann.experiment = 6

        ann.save()
        return HttpResponse(200)
    return HttpResponse(500)


def identify_annotations_results(request, doc_id):
    doc = get_object_or_404(Document, pk=doc_id)

    # If mTurk user not logged in, make a new account for them and set the session
    assignment_id = request.GET.get('assignmentId') #ASSIGNMENT_ID_NOT_AVAILABLE
    hit_id = request.GET.get('hitId')
    # Only available when accepted HIT
    worker_id = request.GET.get('workerId')
    turk_sub_location = request.GET.get('turkSubmitTo')

    return render_to_response('document/concept-recognition-results.jade',
        { 'doc': doc,
          'task_type': 'concept-recognition' },
        context_instance=RequestContext(request))


'''

  Views for completing the Verify Concept task

'''
@login_required
def validate_concepts(request, doc_id):
    doc = get_object_or_404(Document, pk=doc_id)
    relationships = doc.get_conceptrelation_entries_to_validate()

    if len(relationships) is 0:
      return redirect('/document/{0}'.format(doc.pk))

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


    return redirect('/document/'+ str(doc.pk) + '/concepts/' )



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
        return redirect('/document/{0}/results/'.format(doc.pk))

    elif task_type == 'validate-concepts':
        return redirect('/document/{0}/concepts/validate/'.format(doc.pk))
    elif task_type == 'identify-concepts':
        return redirect('/document/{0}/concepts/identify/'.format(doc.pk))

    else:
        doc.update_views(request.user, 'cr', True)
        return redirect('/document/{0}/results/'.format(doc.pk))


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
        return redirect('/document/{0}/'.format(doc.pk))

    elif task_type == 'validate-concepts':
        return redirect('/document/{0}/concepts/validate/'.format(doc.pk))
    elif task_type == 'identify-concepts':
        return redirect('/document/{0}/concepts/identify/'.format(doc.pk))
    else:
        return redirect('/document/{0}/'.format(doc.pk))


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
      return redirect('/document/'+ str(doc.pk) )



class RelationshipTypeViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows users to be viewed or edited.
    """
    queryset = RelationshipType.objects.all()
    serializer_class = RelationshipTypeSerializer


