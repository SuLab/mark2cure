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

# from mark2cure.common.forms import UserForm
from mark2cure.document.models import Document, View, Annotation
from mark2cure.common.forms import MessageForm

from datetime import datetime, timedelta
import math

def home(request):
    if request.user.is_superuser:
      return redirect('/library')

    return render_to_response('landing/index.jade', context_instance=RequestContext(request))

@login_required
def library(request, page_num=1):

    # print "\n\n GET PROFILE"
    # print request.user.profile

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

    '''
      Calc stats for player
    '''
    stats = []
    # Count total seconds spent working
    views = list(View.objects.filter(user = request.user).all().distinct().values_list('updated', 'created', 'section__document'))
    seen = set()
    u_views = [item for item in views if item[2] not in seen and not seen.add(item[2])]
    total_seconds = 0
    for updated, created, document in u_views:
      timediff = (updated - created).total_seconds()
      total_seconds += timediff if (timediff > 3 and timediff < 600) else 0

    stats.append({'t':'Total docs', 'v': len(u_views) })
    stats.append({'t':'Total annotations', 'v': Annotation.objects.filter(view__user = request.user).count()})
    stats.append({'t':'Total time', 'v': str(int(math.floor(total_seconds / 60))) +" mins" })


    return render_to_response('library/index.jade', {
      'docs' : docs,
      'recent': recent_docs,
      'stats': stats}, context_instance=RequestContext(request))


@require_http_methods(["POST"])
def signup(request):
    email = request.POST.get('email', None)
    notify = request.POST.get('email_notify', False)
    if notify is not False:
      notify = True

    if email:
      u, created = User.objects.get_or_create(username=email, email=email)
      if created:
        u.set_password('')
        profile = u.profile
        profile.email_notify = notify
        profile.save()
      u.save()
      return redirect('/')
    return HttpResponse('Unauthorized', status=401)


@require_http_methods(["POST"])
@login_required
def message(request):
    form = MessageForm(request.POST)
    if form.is_valid():
      message = form.save(commit=False)
      message.user = request.user
      message.save()
      return HttpResponse("Success")

    return HttpResponse('Unauthorized', status=401)

def network(request):
    pass


#
# import re
#
# def clean(text):
#   word = re.sub(r'\W+', ' ', text).lower().strip()
#   return word[:-1] if word in ['cancers', 'chordomas'] else word
#
# class Network(Resource):
#   def get(self):
#     #
#     # THIS ALGO IMPLEMENTS ALL USER ANNOTATIONS ACCROSS ALL DOCS
#     #
#     node_list = []
#     link_list = []
#
#     # (TODO) Refine to X # of articles, X days ago, ...
#     annotations = db.session.query(Annotation).filter_by(user = current_user).all()
#
#     # Saving the unique arrays of annotations and documents
#     ann_arr = list(set( [clean(ann.text)  for ann in annotations] ))
#     doc_arr = list(set( [ann.document     for ann in annotations] ))
#
#     # Put the uniq docs at the top of the list then uniq annotation terms after the documents
#     [node_list.append({"name": doc.title, "group": 1})  for doc in doc_arr]
#     [node_list.append({"name": ann, "group": 2})        for ann in ann_arr]
#
#     # For each of the documents with annotations
#     for doc_idx, doc in enumerate( doc_arr ):
#       # Collect all of the cleaned annotations for that document by the current user
#       doc_ann = [clean(ann.text) for ann in doc.annotations if ann.user.id is current_user.id]
#
#       for ann in list(set(doc_ann)):
#         # For this document, connect the doc to the annotation and weight appropriately
#         link_list.append({"source"  : doc_idx,
#                           "target"  : len(doc_arr) + ann_arr.index(ann),
#                           "value"   : doc_ann.count(ann) })
#
#     return jsonify(nodes=node_list, links=link_list)
#
