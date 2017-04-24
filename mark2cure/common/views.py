from django.http import HttpResponse
from django.views.decorators.http import require_http_methods
from django.shortcuts import get_object_or_404, redirect
from django.template.response import TemplateResponse
from django.conf import settings
from django.core.urlresolvers import reverse

from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.messages import get_messages

from allauth.socialaccount.models import SocialApp
from ..userprofile.models import UserProfile
from ..document.models import Annotation, Document, View
from ..task.models import Level, UserQuestRelationship
from ..task.relation.models import Relation

from .utils.mdetect import UAgentInfo
from .forms import SupportMessageForm
from .models import Group

import random


def login_with_zooniverse(request):
    # The view that allows the user to login with their Zooniverse account

    # The public Zooniverse app ID (used for Zoonivser login/authentication) - use the SocialApp DB record to retrieve the app ID
    zooniverseSocialApp = SocialApp.objects.get(provider='zooniverse')
    appId = zooniverseSocialApp.client_id

    ctx = {
            'zooniverse_app_id': appId,
            'zooniverse_callback_url': request.build_absolute_uri(reverse('common:zooniverse-callback'))
            }
    return TemplateResponse(request, 'common/login_with_zooniverse.jade', ctx)

def zooniverse_callback(request):
    # The view that is called after the user logged-in using Zooniverse (this view receives
    # the access token and sends it forward to the allauth module's callback view)
    ctx = {
            'zooniverse_allauth_callback_url': request.build_absolute_uri(reverse('zooniverse_callback'))
            }
    return TemplateResponse(request, 'common/zooniverse_callback.jade', ctx)


def home(request):
    if request.user.is_authenticated():
        return redirect('common:dashboard')

    form = AuthenticationForm()
    quotes = ["To help others.", "Rare disease dad!",
              "In memory of my daughter who had Cystic Fibrosis.",
              "curiosity.", "This is needed.",
              "Goofing off productively.",
              "Community.", "Science!"]
    random.shuffle(quotes)
    groups = Group.objects.all().exclude(name='Practice').order_by('-order')
    ctx = {'form': form,
           'quotes': quotes,
           'groups': groups,
           'ann_count': Annotation.objects.count()}
    return TemplateResponse(request, 'common/landing2.jade', ctx)


def get_started(request):
    uai = UAgentInfo(request.META.get('HTTP_USER_AGENT'), request.META.get('HTTP_ACCEPT'))
    if uai.detectMobileLong():
        return redirect('training:relation-training', part_num=1, step_num=1)

    else:
        if random.random() <= .1:
            return redirect('training:relation-training', part_num=1, step_num=1)
        else:
            return redirect('training:introduction', step_num=1)


def why_mark2cure(request):
    query = UserProfile.objects.exclude(motivation='').order_by('?').values('motivation', 'user')
    return TemplateResponse(request, 'common/why-i-mark2cure.jade', {'profiles': query})


def dashboard(request):
    if not request.user.is_authenticated():
        return redirect('common:home')

    # We redirect user to the training route if their skill is not level 7
    available_tasks = request.user.profile.unlocked_tasks()
    if len(available_tasks) == 0:
        return redirect('training:route')

    welcome = False
    storage = get_messages(request)
    for message in storage:
        if message.message == 'dashboard-unlock-success':
            welcome = True

    uai = UAgentInfo(request.META.get('HTTP_USER_AGENT'), request.META.get('HTTP_ACCEPT'))
    document_ids = list(set(Relation.objects.all().values_list('document', flat=True)))

    er_level = Level.objects.filter(user=request.user, task_type='e').first()
    r_level = Level.objects.filter(user=request.user, task_type='r').first()

    ctx = {'welcome': welcome,
           'mobile': uai.detectMobileLong(),
           'available_tasks': available_tasks,
           'documents': Document.objects.filter(id__in=document_ids)[:100],
           'training_levels': {
               'entity_recognition': er_level.level if er_level else 0,
               'relation': r_level.level if r_level else 0
           },
           'task_stats': {
               'entity_recognition': {
                   'total_score': request.user.profile.score(task='entity_recognition'),
                   'quests_completed': UserQuestRelationship.objects.filter(user=request.user, completed=True).count(),
                   'papers_reviewed': View.objects.filter(user=request.user, completed=True, task_type='cr').count(),
                   'annotations': Annotation.objects.filter(kind='e', view__user=request.user).count()
               },
               'relation': {
                   'total_score': request.user.profile.score(task='relation'),
                   'quests_completed': View.objects.filter(user=request.user, completed=True, task_type='ri').count(),
                   'annotations': Annotation.objects.filter(kind='r', view__user=request.user).count()
               }
           }
           }

    return TemplateResponse(request, 'common/dashboard.jade', ctx)


# removed login required here to allow public to see doc set contributions
def group_view(request, group_stub):
    group = get_object_or_404(Group, stub=group_stub)
    # total annotation count here, plus anns
    top_five, username_list = group.top_five_contributors()

    # (TODO) this should not be hardcoded here; could not open file? Quick Fix. -JF
    # data should come from file in static/data "group_release_dates.txt"
    group_date_dict = {
        "CDG": {"invite": "2015.05.12", "public": "2015.05.21", "closed": "2015.07.29"},
        "alacrima": {"invite": "2015.05.21", "public": "2015.05.22", "closed": "2015.06.19"},
        "OGD": {"invite": "2015.05.29", "public": "2015.05.29", "closed": "2015.11.13"},
        "FBX": {"invite": "2015.06.25", "public": "2015.06.26", "closed": "2015.08.14"},
        "ost": {"invite": "2015.07.31", "public": "2015.08.07", "closed": "2016.04.03"},
        "mfold": {"invite": "2015.09.10", "public": "2015.09.11", "closed": "2015.11.19"},
        "eeyar": {"invite": "2015.11.04", "public": "2015.11.06", "closed": "2015.12.25"},
        "mitomis": {"invite": "2015.11.10", "public": "2015.11.11", "closed": "2016.03.04"},
        "ATGS": {"invite": "2015.12.28", "public": "2015.12.30", "closed": ""},
        "MATG": {"invite": "2016.02.24", "public": "2016.02.26", "closed": ""},
        "MATGS": {"invite": "2016.04.15", "public": "2016.04.15", "closed": ""},
        "training": {"invite": "2015.05.21", "public": "2015.05.21", "closed": ""},
        "HSPT1": {"invite": "2016.04.15", "public": "2016.04.15", "closed": ""},
    }

    try:
        start_date = group_date_dict[group.stub]['invite']
        end_date = group_date_dict[group.stub]['closed']
    except:
        start_date = ""
        end_date = ""

    ctx = {'group': group,
           'top_five': top_five,
           'username_list': username_list,
           'start_date': start_date,
           'end_date': end_date}
    return TemplateResponse(request, 'common/group_home.jade', ctx)


# (TODO) group_network.jade does not exist. Remove this view?
@login_required
def group_network(request, group_stub):
    group = get_object_or_404(Group, stub=group_stub)
    ctx = {'group': group}
    return TemplateResponse(request, 'common/group_network.jade', ctx)


@require_http_methods(['POST'])
def support(request):
    form = SupportMessageForm(data=request.POST)
    if form.is_valid():
        form.save()
        return HttpResponse(200)
    return HttpResponse(500)

