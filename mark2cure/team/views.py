from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User

from django.shortcuts import get_object_or_404, redirect
from django.template.response import TemplateResponse

from mark2cure.userprofile.models import Team

from mark2cure.registration.forms import UserNameChangeForm


def home(request, teamname):
    team = get_object_or_404(Team, name=teamname)
    members = team.userprofile_set.select_related('user')
    ctx = {'team': team, 'members': members}
    return TemplateResponse(request, 'team/public-team.jade', ctx)
