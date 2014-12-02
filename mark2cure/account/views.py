# Create your views here.
from django.shortcuts import render_to_response, get_object_or_404
from django.template import RequestContext
from django.contrib.auth.decorators import login_required
from django.shortcuts import redirect
from django.http import HttpResponse
from django.views.decorators.http import require_http_methods
from django.contrib.auth.models import User
from django.contrib.auth.forms import PasswordChangeForm

from django.contrib.auth import authenticate, login

from mark2cure.account.models import UserProfile
from mark2cure.account.forms import UserForm, UserProfileForm

from rest_framework import viewsets, generics
from rest_framework.response import Response
from rest_framework.decorators import api_view


from mark2cure.common.models import Task, UserQuestRelationship
from brabeion import badges
from brabeion.models import BadgeAward

import datetime
import os


@login_required
def settings(request):
    """"
        Use this endpoint for strict user setting modifications
        like password (auth.user)
    """
    user = request.user


    if request.method == 'POST':
        profileForm = UserProfileForm(request.POST, instance=user.profile)
        if profileForm.is_valid():
            profileForm.save()
            return redirect('mark2cure.account.views.settings')

    profileForm = UserProfileForm(instance=user.profile)
    passwordChangeForm = PasswordChangeForm(user)

    return render_to_response('account/settings.jade',
            {'passwordChangeForm': passwordChangeForm,
             'user_profile': user.profile,
             'user_profile_form': profileForm},
              context_instance=RequestContext(request))


@api_view(['GET'])
@login_required
def user_points(request):
    points_badge = BadgeAward.objects.filter(user=request.user, slug='points').last()
    skill_badge =  BadgeAward.objects.filter(user=request.user, slug='skill').last()

    return Response({
        'points': request.user.userprofile.rating_score,
        'points_level': points_badge.name,
        'skill_level': skill_badge.name
    })


def create(request):
    if request.method == 'POST':
        if request.user.is_authenticated():
            form = UserProfileForm(request.POST, instance=request.user.profile)
            if form.is_valid():
                form.save()
                # They already have an account, let's get them
                # started!
                return redirect('mark2cure.common.views.dashboard')

        else:
            form = UserForm(request.POST)
            if form.is_valid():
                user = form.save()
                profile = user.profile
                profile.save()

                user = authenticate(
                    username=request.POST['username'],
                    password=request.POST['password'])
                login(request, user)

                # In order to create an account they've already done the
                # first 2 training Tasks
                task = Task.objects.first()
                badges.possibly_award_badge("skill_awarded", user=user, level=task.provides_qualification)
                user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
                UserQuestRelationship.objects.create(task=task, user=user, completed=True)

                # Redirect them back b/c of the UserProfileForm
                return redirect('mark2cure.account.views.create')


    if request.user.is_authenticated():
        user_form = UserForm(instance=request.user)
        user_profile_form = UserProfileForm(instance=request.user.profile)
    else:
        user_form = UserForm()
        user_profile_form = UserProfileForm()

    return render_to_response('account/create.jade',
                              {'user_form': user_form,
                               'user_profile_form': user_profile_form},
                              context_instance=RequestContext(request))


@require_http_methods(["POST"])
def newsletter_subscribe(request):
    email = request.POST.get('email', None)
    notify = request.POST.get('email_notify', False)
    if notify is not False:
      notify = True

    if email:
      u, created = User.objects.get_or_create(username=email, email=email)
      if created:
        u.set_password('')
        profile = u.profile
        profile.email_notify = notify
        profile.save()
      u.save()

      return redirect('/')
    return HttpResponse('Unauthorized', status=401)


@login_required
@require_http_methods(['POST'])
def update_profile(request, profile_id):
    """
        User this endpoint for updating profile and game
        specific traits about a profile (auth.user.userprofile)
    """
    profile = get_object_or_404(UserProfile, pk=profile_id)
    form = UserProfileForm(request.POST, instance=profile)

    if form.is_valid():
        profile = form.save()

    return redirect('/account/')
