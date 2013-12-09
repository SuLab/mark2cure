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
from django.contrib.auth import authenticate, login

from mark2cure.document.models import Document, Annotation, View, Section
from mark2cure.document.forms import DocumentForm, AnnotationForm
from mark2cure.document.utils import create_from_pubmed_id
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
      if worker_id:
        return render_to_response('document/read.jade',
                                  { "doc": doc,
                                    "completed": True,
                                    "turk_sub_location": turk_sub_location,
                                    "assignmentId": assignment_id},
                                  context_instance=RequestContext(request))
      else:
        # Move on to another document
        doc = Document.objects.get_random_document()
        return redirect('/document/'+ str(doc.pk) )
    else:

      if worker_id and not request.user.is_authenticated():
        # If it's accepted and a worker that doesn't have an account
        user = get_mturk_account(worker_id)
        user = authenticate(username=user.username, password='')
        login(request, user)

      if request.user.is_authenticated():
        for sec in doc.section_set.all():
          view, created = View.objects.get_or_create(section = sec, user = request.user)

      return render_to_response('document/read.jade',
                                {"doc": doc, "completed": False},
                                context_instance=RequestContext(request))

# showDocument : function(doc_id, assignment_id, hit_id, worker_id, turk_sub) {
#   // console.log('showDocument :: ', doc_id, assignment_id, hit_id, worker_id, turk_sub);
#   User.set('assignment_id', null);
#   window.aws = null;

#   switch(assignment_id) {
#   case undefined:
#     //-- Normal user asking for specific document
#     vent.trigger('navigate:document', {doc_id: doc_id});
#     break;
#   case 'ASSIGNMENT_ID_NOT_AVAILABLE':
#     //-- Preview mode
#     window.aws = {};
#     window.aws.assignment_id = assignment_id;
#     User.set('assignment_id', assignment_id);
#     vent.trigger('navigate:document', {doc_id: doc_id});
#     break;
#   default:
#     window.aws = {};
#     window.aws.turk_sub = turk_sub;
#     window.aws.worker_id = worker_id;
#     window.aws.hit_id = hit_id;
#     window.aws.assignment_id = assignment_id;
#     window.aws.document_id = doc_id;
#     //-- If via AMT, get that user started if not auth'd already
#     User.set('assignment_id', assignment_id);

#     if( User.authenticated() && User.get('mturk') ) {
#       vent.trigger('navigate:document', {doc_id: doc_id});
#     } else {
#       User.set('username', worker_id);
#       User.set('mturk', true);
#       User.save(null, {success: function() {
#         //-- After our user is saved, go ahead to get the document
#         vent.trigger('navigate:document', {doc_id: doc_id});
#       }});
#     }
#     break;
#   }

# },

# views = db.session.query(View).filter_by(user = current_user).filter_by( document = document ).all()
# if len(views):
  # if current_user.mturk:
    # t = Turk()
    # t.mtc.block_worker(current_user.username, "Attempted to submit same document multiple times.")
  # raise ValueError("Cannot submit a document twice")

# if document.validate and current_user.mturk:
#       # If this is a validate document, check the user's history, if it's their 3rd submission
#       # or more run test to potentially fail if poor performance
#       valid_views = db.session.query(View).\
#           filter_by(user = current_user).\
#           filter( View.document.has(validate=1) ).\
#           order_by( desc(View.created) ).\
#           limit(3).all()
#       if len(valid_views) is 3:
#         if sum(1 for x in valid_views if gold_matches(x.user, x.document) >= 1) is not 3:
#           print "failed"
#           t = Turk()
#           t.mtc.block_worker(current_user.username, "Failed to properly answer golden master performance documents")
#


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

      ann.save()
      return HttpResponse("Success")

    return HttpResponse("Failed")

