'''
Use these views to serve up static pages,
e.g. an about, FAQ, or help page.
'''

from django.contrib.auth.decorators import login_required
from django.template import RequestContext
from django.shortcuts import get_object_or_404, render_to_response, redirect
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.views.decorators.http import require_http_methods
from django.contrib.auth.models import User
from django.shortcuts import redirect
from django.http import HttpResponse

from mark2cure.common.forms import UserForm
from mark2cure.document.models import Document

from datetime import datetime, timedelta

def home(request):
    if request.user.is_superuser:
      return redirect('/library')

    return render_to_response('landing/index.jade', context_instance=RequestContext(request))

@login_required
def library(request, page_num=1):

    query = request.GET.get('search')
    if query:
      doc_list = Document.objects.filter(section__text__search = query).distinct()
    else:
      doc_list = Document.objects.all()

    doc_list_paginator = Paginator(doc_list, 18)
    try:
        docs = doc_list_paginator.page(page_num)
    except PageNotAnInteger:
        docs = doc_list_paginator.page(1)
    except EmptyPage:
        docs = doc_list_paginator.page(paginator.num_pages)

    recent_docs = Document.objects\
        .filter(section__view__user = request.user)\
        .order_by('-created')\
        .distinct()[3:]

    return render_to_response('library/index.jade', {'docs' : docs, 'recent': recent_docs}, context_instance=RequestContext(request))


@require_http_methods(["POST"])
def create(request):
    form = UserForm(request.POST)
    if form.is_valid():
      user = form.save()
      return redirect('/')
    else:
      form = UserForm()
      return HttpResponse('Unauthorized', status=401)
