from django.contrib.auth.decorators import login_required
from django.views.decorators.cache import cache_page
from django.contrib.auth.models import User

from django.shortcuts import get_object_or_404, redirect
from django.template.response import TemplateResponse

from allauth.account.forms import UserForm

from .forms import UserProfileForm, TeamForm
from ..task.models import Level

from rest_framework.decorators import api_view
from rest_framework.response import Response
from django.contrib import messages


@cache_page(60 * 15)
def public_profile(request, username):
    user = get_object_or_404(User, username=username)

    ctx = {'player': user,
           'owner': True if request.user == user else False,
           'contributed_groups': user.profile.contributed_groups()}

    return TemplateResponse(request, 'userprofile/public-profile.html', ctx)


@login_required
def settings(request):
    '''
        Use this endpoint for strict user setting modifications
        like password (auth.user)
    '''
    user_change_form = UserForm(user=request.user, data=request.POST or None)

    user_profile_form = UserProfileForm(instance=request.user.profile, data=request.POST or None)
    team_form = TeamForm(instance=None, data=request.POST or None)

    if request.method == 'POST':
        UserForm(user=request.user, data=request.POST or None)
        profile = user_profile_form.save()

        if team_form.is_valid() and team_form.instance.name != '':
            team = team_form.save(commit=False)
            team.owner = request.user
            team.save()

            profile.team = team
            profile.save()

        messages.info(request, '<p class="lead text-center my-1">Profile Successfully Updated</p>', extra_tags='safe alert-info')
        return redirect('userprofile:settings')

    ctx = {'user_change_form': user_change_form,
           'team_form': team_form,
           'user_profile_form': user_profile_form}
    return TemplateResponse(request, 'userprofile/settings.html', ctx)


@api_view(['GET'])
@login_required
def user_points(request):

    return Response({
        'points': request.user.profile.score(),
        # (TODO) implement using new score method
        'points_level': 'Hard Worker',
    })
