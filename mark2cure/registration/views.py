from django.conf import settings
from django.template.response import TemplateResponse
from django.shortcuts import redirect

from django.contrib.auth.decorators import login_required
from django.contrib.auth import authenticate, login as auth_login, logout as auth_logout
from django.contrib import messages

from ..task.models import Level, Task, UserQuestRelationship
from mark2cure.userprofile.forms import UserProfileForm
from mark2cure.score.models import Point

from django.core.urlresolvers import reverse
from django.contrib.auth.views import password_reset, password_reset_confirm
from django.utils.translation import ugettext as _

from django.utils import timezone
from datetime import datetime
from urllib import urlencode
from brabeion import badges
import forms
import utils

import os


def user_creation(request):
    if request.user.is_authenticated():
        return redirect('common:dashboard')

    user_create_form = forms.UserCreateForm(data=request.POST or None)
    if request.POST and user_create_form.is_valid():
        user = user_create_form.save()

        user = authenticate(
            username=request.POST['username'],
            password=request.POST['password1'])
        auth_login(request, user)

        # (TODO) Needs to know what training they're coming from
        task = Task.objects.first()
        Level.objects.create(user=user, task_type='e', level=task.provides_qualification, created=timezone.now())

        from django.contrib.contenttypes.models import ContentType
        content_type = ContentType.objects.get_for_model(task)
        Point.objects.create(user=request.user, amount=task.points, content_type=content_type, object_id=task.id, created=timezone.now())

        UserQuestRelationship.objects.create(task=task, user=user, completed=True)

        # Redirect them back b/c of the UserProfileForm
        return redirect('registration:user_creation_settings')

    return TemplateResponse(request, 'registration/create.jade', {'form': user_create_form})


@login_required
def user_creation_settings(request):
    user_change_form = forms.UserNameChangeForm(instance=request.user, data=request.POST or None)
    user_profile_form = UserProfileForm(instance=request.user.profile, data=request.POST or None)

    if request.POST and user_profile_form.is_valid():
        user_change_form.save()
        user_profile_form.save()
        # They already have an account, let's get them
        # started!
        return redirect('common:dashboard')

    return TemplateResponse(request, 'registration/create-settings.jade',
            {'user_change_form': user_change_form,
             'user_profile_form': user_profile_form})


def reset_confirm(request, uidb36=None, token=None):
    return password_reset_confirm(request,
            template_name='password-reset/password_reset_confirm.jade',
            uidb36=uidb36,
            token=token,
            post_reset_redirect=reverse('common:dashboard'))


def reset(request):
    return password_reset(request,
            template_name='password-reset/password_reset_complete.jade',
            email_template_name='app/reset_email.html',
            subject_template_name='app/reset_subject.txt',
            post_reset_redirect=reverse('app:login'))
