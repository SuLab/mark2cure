from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User

from django.shortcuts import get_object_or_404, redirect
from django.template.response import TemplateResponse

from .forms import UserProfileForm
from mark2cure.registration.forms import UserNameChangeForm
from rest_framework.decorators import api_view
from rest_framework.response import Response
from django.contrib import messages

from brabeion.models import BadgeAward

from mark2cure.document.models import Annotation
from mark2cure.common.models import UserQuestRelationship


def public_profile(request, username):
    user = get_object_or_404(User, username=username)
    anns = Annotation.objects.filter(view__user=user)
    quests = UserQuestRelationship.objects.filter(user=user, completed=True)
    # (TODO) Latest Qoutes

    ctx = {'player': user,
           'owner': True if request.user == user else False,
           'quests_count': quests.count(),
           'annotations_count': anns.count()}

    return TemplateResponse(request, 'userprofile/public-profile.jade', ctx)

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
        messages.info(request, 'helloooo', extra_tags='safe alert-success')
        return redirect('userprofile:settings')

    ctx = {'user_change_form': user_change_form,
           'user_profile_form': user_profile_form}
    return TemplateResponse(request, 'userprofile/settings.jade', ctx)


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


