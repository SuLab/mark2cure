'''
Use these views to serve up static pages,
e.g. an about, FAQ, or help page.
'''

from django.contrib.auth.decorators import login_required
from django.template import RequestContext
from django.shortcuts import get_object_or_404, render_to_response, redirect
from django.core.paginator import Paginator
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
def library(request):
    doc_list = Document.objects.all()
    docs = Paginator(doc_list, 20).page(1)
    return render_to_response('library/index.jade', {'docs' : docs}, context_instance=RequestContext(request))

@require_http_methods(["POST"])
def create(request):
    print request.POST
    form = UserForm(request.POST)
    if form.is_valid():
      user = form.save()
      return redirect('/')
    else:
      form = UserForm()
      return HttpResponse('Unauthorized', status=401)
