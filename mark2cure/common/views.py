from django.http import HttpResponse
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms import AuthenticationForm
from django.shortcuts import get_object_or_404, redirect
from django.template.response import TemplateResponse
from django.views.decorators.cache import cache_page
from django.contrib.messages import get_messages
from django.contrib import messages

from ..userprofile.models import UserProfile
from .models import Group
from .forms import SupportMessageForm

import random


#@cache_page(62 * 60)
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
    return TemplateResponse(request, 'common/why-mark2cure.jade', {'profiles': query})


@login_required
def dashboard(request):
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


@login_required
def group_view(request, group_stub):
    group = get_object_or_404(Group, stub=group_stub)
    ctx = {'group': group}
    return TemplateResponse(request, 'common/group_home.jade', ctx)


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


