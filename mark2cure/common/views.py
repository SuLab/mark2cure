from django.http import HttpResponse
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms import AuthenticationForm
from django.shortcuts import get_object_or_404, redirect
from django.template.response import TemplateResponse
from django.contrib.messages import get_messages
from django.contrib import messages

from ..userprofile.models import UserProfile
from .models import Group
from .forms import SupportMessageForm

import random


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
    ctx = {'form': form, 'quotes': quotes, 'groups': groups}
    return TemplateResponse(request, 'common/landing2.jade', ctx)


def beta(request):
    return redirect('common:home')


def why_mark2cure(request):
    query = UserProfile.objects.exclude(motivation='').order_by('?').values('motivation', 'user')
    return TemplateResponse(request, 'common/why-i-mark2cure.jade', {'profiles': query})


@login_required
def dashboard(request):
    # We redirect user to the training route if their skill is not level 7
    if not request.user.profile.highest_level("skill").level == 7:
        return redirect('training:route')

    welcome = False
    storage = get_messages(request)
    for message in storage:
        if message.message == 'dashboard-unlock-success':
            welcome = True

    msg = '<p class="lead text-center">Click on one of the quest numbers below to start the quest. Your contributions are important so complete as many quests as you can.</p>'
    messages.info(request, msg, extra_tags='safe alert-success')

    ctx = {'welcome': welcome}
    return TemplateResponse(request, 'common/dashboard.jade', ctx)

# removed login required here to allow public to see doc set contributions
def group_view(request, group_stub):
    group = get_object_or_404(Group, stub=group_stub)
    # total annotation count here, plus anns
    top_five, username_list = group.top_five_contributors()

    # (TODO) this should not be hardcoded here; could not open file? Quick Fix. -JF
    # data should come from file in static/data "group_release_dates.txt"
    group_date_dict = {
    "CDG": {"invite":"2015.05.12", "public": "2015.05.21", "closed":"2015.07.29"},
    "alacrima": {"invite":"2015.05.21", "public":"2015.05.22", "closed":"2015.06.19"},
    "OGD": {"invite":"2015.05.29","public":"2015.05.29", "closed":"2015.11.13"},
    "FBX": {"invite":"2015.06.25","public":"2015.06.26","closed":"2015.08.14"},
    "ost": {"invite":"2015.07.31",
    "public": "2015.08.07","closed": "2016.04.03"},
    "mfold": {"invite":"2015.09.10", "public": "2015.09.11","closed": "2015.11.19"},
    "eeyar": {"invite": "2015.11.04","public":"2015.11.06", "closed":"2015.12.25",},
    "mitomis": {"invite":"2015.11.10", "public": "2015.11.11", "closed": "2016.03.04"},
    "ATGS": {"invite":"2015.12.28", "public": "2015.12.30", "closed": ""},
    "MATG": {"invite":"2016.02.24", "public": "2016.02.26", "closed": ""},
    "MATGS": {"invite": "", "public": "", "closed": ""},
    "training": {"invite": "", "public": "", "closed": ""},
    "HSPT1": {"invite": "", "public": "", "closed": ""},
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


#(TODO) group_network.jade does not exist. Remove this view?
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

