from django.template.response import TemplateResponse
from django.shortcuts import get_object_or_404

from ..userprofile.models import Team


def home(request, slug):
    team = get_object_or_404(Team, slug=slug)
    members = team.userprofile_set.select_related('user')
    ctx = {'team': team, 'members': members}
    return TemplateResponse(request, 'team/public-team.html', ctx)
