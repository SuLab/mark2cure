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
    return TemplateResponse(request, 'common/landing2.jade', {'form': form, 'quotes': quotes})


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
    ctx = {'group': group,
           'top_five': top_five,
           'username_list': username_list}
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

