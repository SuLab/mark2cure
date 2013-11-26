'''
Use these views to serve up static pages,
e.g. an about, FAQ, or help page.
'''

from django.contrib.auth.decorators import login_required
from django.template import RequestContext
from django.shortcuts import get_object_or_404, redirect
from django.shortcuts import render_to_response
from django.core.paginator import Paginator

from mark2cure.document.models import Document

from datetime import datetime, timedelta

# @login_required
def home(request):
    # doc_list = Document.objects.all()
    # docs = Paginator(doc_list, 20).page(1)

    return render_to_response('landing/index.jade', context_instance=RequestContext(request))
