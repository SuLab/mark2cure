from django.contrib.auth.decorators import login_required
from django.template import RequestContext
from django.shortcuts import get_object_or_404, render_to_response, redirect
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.views.decorators.http import require_http_methods
from django.contrib.auth.models import User
from django.shortcuts import redirect
from django.http import HttpResponse

from mark2cure.document.models import Document, View, Annotation
from mark2cure.common.forms import MessageForm
from mark2cure.common.models import SurveyFeedback

from datetime import datetime, timedelta
import math


def home(request):
    if request.user.is_superuser:
      return redirect('/library')

    return render_to_response('landing/index.jade', context_instance=RequestContext(request))


def mturk(request):
    '''
      All AWS people will come in through here
      document app will *always* require authentication
    '''
    assignment_id = request.GET.get('assignmentId') #ASSIGNMENT_ID_NOT_AVAILABLE
    worker_id = request.GET.get('workerId')
    turk_sub_location = request.GET.get('turkSubmitTo')

    # If mTurk user not logged in, make a new account for them and set the session
    if assignment_id == 'ASSIGNMENT_ID_NOT_AVAILABLE':
        logout(request)

        doc = Document.objects.get(pk=278)
        return render_to_response('document/concept-recognition.jade',
                                  { 'doc': doc,
                                    'task_type': 'concept-recognition' },
                                  context_instance=RequestContext(request))

    # If they've accepted a HIT
    if worker_id and not request.user.is_authenticated():
        # If it's accepted and a worker that doesn't have an account
        # Make one and log them in
        user = get_mturk_account(worker_id)
        user = authenticate(username=user.username, password='')
        login(request, user)


    if request.user.is_authenticated():
        user_profile = request.user.userprofile

        if assignment_id and turk_sub_location and worker_id:
            user_profile.turk_submit_to = turk_sub_location
            user_profile.turk_last_assignment_id = assignment_id
            user_profile.save()

        if user_profile.softblock:
            return redirect('mark2cure.common.views.softblock')


    # (TODO) Some magic document selection here
    doc_id = 278
    return redirect('mark2cure.document.views.identify_annotations', doc_id)



def softblock(request):
    return render_to_response('common/softblock.jade', context_instance=RequestContext(request))


@login_required
def library(request, page_num=1):

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
        .distinct()[:3]

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


@require_http_methods(["POST"])
@login_required
def survey(request):
    for k, v in request.POST.iteritems():
        print k, v
        if(k != "csrfmiddlewaretoken"):
          sf = SurveyFeedback(question = k, response = v, user = request.user)
          sf.save()

    return HttpResponse("Success")
    #return HttpResponse('Unauthorized', status=401)


