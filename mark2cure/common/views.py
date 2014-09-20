from django.template import RequestContext
from django.shortcuts import get_object_or_404, render_to_response, redirect
from django.views.decorators.http import require_http_methods

from django.contrib.auth.models import User
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.decorators import login_required
from django.contrib.auth import authenticate, login, logout

from django.http import HttpResponse
from django.conf import settings
from django.core.mail import send_mail
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

from mark2cure.document.models import Document, View, Annotation, Activity
from zinnia.models import Entry

from datetime import datetime, timedelta
import math, random, logging
logger = logging.getLogger(__name__)


def home(request):

    if request.user.is_authenticated():
      return redirect('mark2cure.common.views.dashboard')

    form = AuthenticationForm()
    return render_to_response('common/index.jade',
                              {'form': form},
                              context_instance=RequestContext(request))


def introduction(request):
    return render_to_response('training/basics.jade',
                              {}, context_instance=RequestContext(request))


def training_read(request):
    return render_to_response('training/training_read.jade',
                              {}, context_instance=RequestContext(request))

def training_one(request, step_num):
    if step_num == 'complete':
        return render_to_response(
            'training/intro-1/complete.jade',
            {}, context_instance=RequestContext(request))

    step_num = int(step_num)
    next_ = step_num+1
    if step_num == 1:
        header1 = "Let's start by marking diseases"
        header2 = "Mark all the disease terms in the sentence below."
        paragraph = "Does choice of insulin regimen really matter in the management of diabetes?"
        answers = ['diabetes']

    if step_num == 2:
        header1 = "Sometimes you will see multiple instances of the same disease - Be sure to mark them all!"
        header2 = "Mark all the disease terms in the sentence below."
        paragraph = "To assess the management of diabetes, we reviewed records of 20 diabetes patients."
        answers = ['diabetes', 'diabetes']

    if step_num == 3:
        header1 = "Sometimes the disease is described by a conjuction of several words. Mark these disease conjunctions as a single span of text."
        header2 = "Mark all the disease terms in the sentence below."
        answers = ['diabetes', 'type 2 diabetes mellitus'];
        paragraph = "Of the 20 diabetes patients, 17 had type 2 diabetes mellitus."

    if step_num == 4:
        header1 = "Sometimes the disease conjunctions are separated by words like 'and/or'. Decide if there are two distinct diseases to highlight, or if it is a single disease conjunction."
        header2 = "Mark all the disease terms in the sentence below."
        paragraph = "The remaining 3 had inherited and/or type I diabetes mellitus."
        answers = ['inherited and/or type I diabetes mellitus'];

    if step_num == 5:
        header1 = "Sometimes different diseases are discussed. Mark all the diseases below."
        header2 = "Remember to mark disease conjunction as spans and to mark different diseases separately."
        paragraph = "Of the 20 patients, 10 patients were also diagnosed with heart disease or rheumatoid arthritis."
        answers = ['heart disease', 'rheumatoid arthritis'];

    if step_num == 6:
        header1 = "Sometimes the disease terms are abbreviated. Mark all instances of disease abbreviations."
        header2 = "Mark all disease terms in the sentence below."
        paragraph = "We will discuss the effect of different insulin regimen on type 2 diabetes mellitus patients (ie- T2DM patients) ..."
        answers = ['type 2 diabetes mellitus', 'T2DM'];

    if step_num == 7:
        header1 = "Practice what you've learned so far..."
        header2 = "Try marking the disease and disease abbreviations in this phrase now!"
        paragraph = "... with or without rheumatoid arthritis ( RA ) or heart disease ( HD )."
        answers = ['rheumatoid arthritis', 'RA', 'heart disease', 'HD'];

    if step_num == 8:
        header1 = "Now let's mark some symptoms."
        header2 = "Symptoms are the physical manifestations of the disease. Mark the symptoms in the sentence below."
        paragraph = "In particular, we will examine the effects of these insulin regimen on the symptoms of the diseases. We will focus on some common symptoms such as fatigue as well as ..."
        answers = ['fatigue'];

    if step_num == 9:
        header1 = "Sometimes a single symptom is described using more than one word."
        header2 = "Mark these symptoms as single spans of text. Mark the symptom in the text below."
        paragraph = "... as well as frequent urination. Another crucial symptom ..."
        answers = ['frequent urination'];

    if step_num == 10:
        header1 = "Sometimes a single symptom is described with a long block of text and may have joining terms such as 'and' or 'or'."
        header2 = "In that case, highlihgt the entire symptom as a single span of text. Finish Training #1 by marking the symptom in the text below."
        paragraph = "Another crucial symptom that we will investigate includes tingling sensations and/or numbness in the hands or feet."
        answers = ['tingling sensations and/or numbness in the hands or feet'];
        next_ = 'complete'

    return render_to_response(
        'training/intro-1/read.jade',
        {'training_num': 1,
         'step_num': step_num,
         'header1': header1,
         'header2': header2,
         'paragraph': paragraph,
         'answers': answers,
         'next': next_},
        context_instance=RequestContext(request))

@login_required
def training_two(request, step_num):
    if step_num == 'complete':
        profile = request.user.profile
        profile.training_complete = True
        profile.save()

        return redirect('mark2cure.common.views.dashboard')

    return render_to_response(
        'training/intro-2/step-{step_num}.jade'.format(step_num=step_num),
        {'step_num': step_num},
        context_instance=RequestContext(request))


@login_required
def dashboard(request):
    if not request.user.profile.training_complete:
        return redirect('mark2cure.common.views.training_read')

    profile = request.user.profile
    posts = Entry.objects.all()[:3]
    return render_to_response('common/dashboard.jade',
                              {'posts': posts,
                               'profile': profile},
                              context_instance=RequestContext(request))


@login_required
def profile_survey(request):
    user_profile = request.user.userprofile

    if request.method == 'POST':
        form = ProfileSurveyForm(request.POST, instance = user_profile)
        if form.is_valid():
            form.save()
            return redirect('mark2cure.common.views.library')

    else:
        form = ProfileSurveyForm(instance = user_profile)

    return render_to_response('common/survey.jade',
                              {'user_profile': user_profile, 'form': form },
                              context_instance=RequestContext(request))


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
        # print k, v
        if(k != "csrfmiddlewaretoken"):
          sf = SurveyFeedback(question = k, response = v, user = request.user)
          sf.save()

    return HttpResponse("Success")


