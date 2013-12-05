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

from mark2cure.document.models import Document, Annotation, View, Section
from mark2cure.document.forms import DocumentForm, AnnotationForm
from mark2cure.document.utils import create_from_pubmed_id
from mark2cure.common.utils import get_timezone_offset

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


@login_required
def read(request, doc_id):
    if request.method == 'POST':
      # Move on
      doc = Document.objects.get_random_document()
      return redirect('/document/'+ str(doc.pk) )
    else:
      doc = get_object_or_404(Document, pk=doc_id)
      for sec in doc.section_set.all():
        view, created = View.objects.get_or_create(section = sec, user = request.user)

      return render_to_response('document/read.jade',
                                {"doc": doc},
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

      ann.save()

      return HttpResponse("Success")
    return HttpResponse("Failed")

