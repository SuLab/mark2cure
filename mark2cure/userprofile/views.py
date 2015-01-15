from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.contrib.auth.forms import PasswordChangeForm

from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect, Http404
from django.shortcuts import get_object_or_404, redirect
from django.template.response import TemplateResponse
from django.views.decorators.http import require_POST
from django.utils.translation import ugettext as _

from .forms import UserProfileForm
from mark2cure.registration.forms import UserNameChangeForm
from rest_framework.decorators import api_view


@login_required
def settings(request):
    '''
        Use this endpoint for strict user setting modifications
        like password (auth.user)
    '''
    user_change_form = UserNameChangeForm(instance=request.user, data=request.POST or None)
    user_profile_form = UserProfileForm(instance=request.user.profile, data=request.POST or None)

    if request.method == 'POST':
        user_change_form.save()
        user_profile_form.save()
        return redirect('account:settings')

    ctx = { 'user_change_form': user_change_form,
            'user_profile_form': user_profile_form}
    return TemplateResponse(request, 'account/settings.jade', ctx)


@api_view(['GET'])
@login_required
def user_points(request):
    points_badge = BadgeAward.objects.filter(user=request.user, slug='points').last()
    skill_badge = BadgeAward.objects.filter(user=request.user, slug='skill').last()

    return Response({
        'points': request.user.userprofile.rating_score,
        'points_level': points_badge.name,
        'skill_level': skill_badge.name
    })


