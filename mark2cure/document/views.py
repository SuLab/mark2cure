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

from mark2cure.document.models import Document, Annotation, View, Section
from mark2cure.document.forms import DocumentForm, AnnotationForm
from mark2cure.document.utils import create_from_pubmed_id, check_validation_status
from mark2cure.common.utils import get_timezone_offset, get_mturk_account

from copy import copy

import oauth2 as oauth
import json


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




# @login_required
def read(request, doc_id):
    # If they're attempting to view or work on the document
    doc = get_object_or_404(Document, pk=doc_id)

    # If mTurk user not logged in, make a new account for them and set the session
    assignment_id = request.GET.get('assignmentId') #ASSIGNMENT_ID_NOT_AVAILABLE
    hit_id = request.GET.get('hitId')
    # Only available when accepted HIT
    worker_id = request.GET.get('workerId')
    turk_sub_location = request.GET.get('turkSubmitTo')

    if request.method == 'POST':
      '''
        If the user if submitting results for a document an document and sections
      '''

      if worker_id:

        if request.user.is_authenticated():
          check_validation_status(request.user, doc)

        return render_to_response('document/read.jade',
                                  { "doc": doc,
                                    "completed": True,
                                    "turk_sub_location": turk_sub_location,
                                    "instruct_bool": "block" if assignment_id == "ASSIGNMENT_ID_NOT_AVAILABLE" else "none",
                                    "assignmentId": assignment_id},
                                  context_instance=RequestContext(request))
      else:

        if request.user.is_authenticated():
          # Update the timestamps
          for sec in doc.section_set.all():
            view = get_object_or_404(View, section = sec, user = request.user)
            view.save()

        # Move on to another document
        doc = Document.objects.get_random_document()
        return redirect('/document/'+ str(doc.pk) )

    else:
      '''
        If the user is fetching / viewing a document and sections
      '''
      completed = False

      if assignment_id == "ASSIGNMENT_ID_NOT_AVAILABLE":
        logout(request)

      if worker_id and not request.user.is_authenticated():
        # If it's accepted and a worker that doesn't have an account
        user = get_mturk_account(worker_id)
        user = authenticate(username=user.username, password='')
        login(request, user)

      if request.user.is_authenticated():
        sections = doc.section_set.all()
        for sec in sections:
          view, created = View.objects.get_or_create(section = sec, user = request.user)

        gen_u_view = View.objects.filter(section__document = doc, user = request.user).values('updated', 'created').first()
        timediff = (gen_u_view['updated'] - gen_u_view['created']).total_seconds()
        completed = timediff > 2

        if completed:
          results = {}
          for sec in sections:
            print " ~ ~ ~ ~ ~ ~ ~ "
            anns = Annotation.objects.filter(view__section = sec).values_list('text', 'start').all()
            print anns

      return render_to_response('document/read.jade',
                                { "doc": doc,
                                  "completed": False,
                                  "instruct_bool": "block" if assignment_id == "ASSIGNMENT_ID_NOT_AVAILABLE" else "none",
                                  "assignmentId": assignment_id},
                                context_instance=RequestContext(request))



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


@login_required
@require_http_methods(["POST"])
def createannotation(request, doc_id, section_id):
    section = get_object_or_404(Section, pk=section_id)
    view, created = View.objects.get_or_create(section = section, user = request.user)

    form = AnnotationForm(request.POST)
    if form.is_valid():
        ann = form.save(commit=False)

        ann.view = view
        ann.type = "disease"
        ann.user_agent = request.META['HTTP_USER_AGENT']
        ann.player_ip = request.META['REMOTE_ADDR']

        if request.user.profile.mturk:
          ann.experiment = 3

        ann.save()
        return HttpResponse(200)
    return HttpResponse(500)

