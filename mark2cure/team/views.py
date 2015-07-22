from django.template.response import TemplateResponse
from django.shortcuts import get_object_or_404

from ..userprofile.models import Team


def home(request, teamname):
    team = get_object_or_404(Team, name=teamname)
    members = team.userprofile_set.select_related('user')
    ctx = {'team': team, 'members': members}
    return TemplateResponse(request, 'team/public-team.jade', ctx)
