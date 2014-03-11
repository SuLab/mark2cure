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
from mark2cure.document.utils import update_views, create_from_pubmed_id, check_validation_status
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
                              {"docs": docs},
                              context_instance=RequestContext(request))




def identify_annotations(request, doc_id):
    # If they're attempting to view or work on the document
    doc = get_object_or_404(Document, pk=doc_id)

    # If mTurk user not logged in, make a new account for them and set the session
    assignment_id = request.GET.get('assignmentId') #ASSIGNMENT_ID_NOT_AVAILABLE
    hit_id = request.GET.get('hitId')
    # Only available when accepted HIT
    worker_id = request.GET.get('workerId')
    turk_sub_location = request.GET.get('turkSubmitTo')

    '''
      If the user is fetching / viewing a document and sections
    '''
    completed = False

    if len(doc.section_set.filter(kind="a")) is 0:
      return redirect('/document/'+ str(doc.pk) + '/concepts/validate/' )


    if assignment_id == "ASSIGNMENT_ID_NOT_AVAILABLE":
      logout(request)

    if worker_id and not request.user.is_authenticated():
      # If it's accepted and a worker that doesn't have an account
      user = get_mturk_account(worker_id)
      user = authenticate(username=user.username, password='')
      login(request, user)

    if request.user.is_authenticated():
      update_views(request.user, doc, 'cr')

      gen_u_view = View.objects.filter(task_type = 'cr', section__document = doc, user = request.user).values('updated', 'created').first()
      timediff = (gen_u_view['updated'] - gen_u_view['created']).total_seconds()
      completed = timediff > 2

    return render_to_response('document/concept-recognition.jade',
                              { "doc": doc,
                                "completed": False,
                                "instruct_bool": "block" if assignment_id == "ASSIGNMENT_ID_NOT_AVAILABLE" else "none",
                                "assignmentId": assignment_id},
                              context_instance=RequestContext(request))

@login_required
@require_http_methods(["POST"])
def create_annotation(request, doc_id, section_id):
    '''
      This is broken out because there can be many submissions per document
      We don't want to use these submission to direct the user to elsewhere in the app
    '''
    section = get_object_or_404(Section, pk=section_id)
    view, created = View.objects.get_or_create(task_type = "cr", section = section, user = request.user)

    form = AnnotationForm(request.POST, view)
    if form.is_valid():
        ann = form.save(commit=False)

        ann.view = view
        ann.type = "disease"
        ann.user_agent = request.META['HTTP_USER_AGENT']
        ann.player_ip = request.META['REMOTE_ADDR']

        if request.user.profile.mturk:
          ann.experiment = 5

        ann.save()
        return HttpResponse(200)
    return HttpResponse(500)



@login_required
def validate_concepts(request, doc_id):
    doc = get_object_or_404(Document, pk=doc_id)
    relationships = doc.get_conceptrelation_entries_to_validate()

    if len(relationships) is 0:
      return redirect('/document/{0}'.format(doc.pk))

    return render_to_response('document/verify-relationships.jade',
        {'doc': doc, 'relationships' : relationships },
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



@login_required
def identify_concepts(request, doc_id):
    # If they're attempting to view or work on the document
    doc = get_object_or_404(Document, pk=doc_id)
    concepts = doc.get_concepts_for_classification()
    return render_to_response('document/identify-concepts.jade',
        {'doc': doc, 'concepts': concepts[:3]},
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


@require_http_methods(['POST'])
def next(request, doc_id):
    '''
      If the user if submitting results for a document an document and sections
    '''
    doc = get_object_or_404(Document, pk=doc_id)

    # # If mTurk user not logged in, make a new account for them and set the session
    # assignment_id = request.GET.get('assignmentId') #ASSIGNMENT_ID_NOT_AVAILABLE
    # hit_id = request.GET.get('hitId')
    # # Only available when accepted HIT
    # worker_id = request.GET.get('workerId')
    # turk_sub_location = request.GET.get('turkSubmitTo')

    # if worker_id:

    #   if request.user.is_authenticated():
    #     check_validation_status(request.user, doc)
    #     # Update the timestamps
    #     update_views(request.user, doc)
    #   return render_to_response('document/concept-recognition.jade',
    #                             { "doc": doc,
    #                               "completed": True,
    #                               "turk_sub_location": turk_sub_location,
    #                               "instruct_bool": "block" if assignment_id == "ASSIGNMENT_ID_NOT_AVAILABLE" else "none",
    #                               "assignmentId": assignment_id},
    #                             context_instance=RequestContext(request))
    # else:

    if request.user.is_authenticated():
      # Update the timestamps
      update_views(request.user, doc)

    kind = request.POST.get('task_type')
    # Move on to another document
    doc = Document.objects.get_random_document()

    if kind == 'concept-recognition':
        return redirect('/document/{0}/'.format(doc.pk))
    elif kind == 'validate-concepts':
        return redirect('/document/{0}/concepts/validate/'.format(doc.pk))
    elif kind == 'identify-concepts':
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


def annotation_results(request, doc_id):
    doc = get_object_or_404(Document, pk=doc_id)



    anns = Annotation.objects.filter(
        view__section__document = doc,
        kind = 'e')

    return render_to_response('document/concept-recognition-results.jade',
        { "doc": doc },
                              context_instance=RequestContext(request))


class RelationshipTypeViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows users to be viewed or edited.
    """
    queryset = RelationshipType.objects.all()
    serializer_class = RelationshipTypeSerializer


