from django.conf import settings
from django.template.response import TemplateResponse
from django.shortcuts import redirect

from django.contrib.auth.decorators import login_required
from django.contrib.auth import authenticate, login as auth_login
from django.contrib import messages

from mark2cure.common.models import Task, UserQuestRelationship
from mark2cure.userprofile.forms import UserProfileForm

from django.core.urlresolvers import reverse
from django.contrib.auth.views import password_reset, password_reset_confirm
from django.utils.translation import ugettext as _

from brabeion import badges
from .models import EmailConfirmationRequest, EmailChangeRequest
from urllib import urlencode

import forms
import utils

import os


def user_creation(request):
    user_create_form = forms.UserCreateForm(data=request.POST or None)
    if request.POST and user_create_form.is_valid():
        user = user_create_form.save()

        user = authenticate(
            username=request.POST['username'],
            password=request.POST['password1'])
        auth_login(request, user)

        # In order to create an account they've already done the
        # first 2 training Tasks
        task = Task.objects.first()
        badges.possibly_award_badge('skill_awarded', user=user, level=task.provides_qualification)
        user.profile.rating.add(score=task.points, user=None, ip_address=os.urandom(7).encode('hex'))
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


def request_email_confirmation(request):
    local_host = utils.get_local_host(request)
    form = forms.RequestEmailConfirmationForm(local_host=local_host,
                                              data=request.POST or None)
    if form.is_valid():
        form.send()
        msg = _('Confirmation email has been sent. '
                'Please check your inbox.')
        messages.success(request, msg)
        return redirect(settings.LOGIN_REDIRECT_URL)

    return TemplateResponse(request,
                            'registration/request_email_confirmation.jade',
                            {'form': form})


@login_required
def request_email_change(request):
    form = forms.RequestEmailChangeForm(
        local_host=utils.get_local_host(request), user=request.user,
        data=request.POST or None)
    if form.is_valid():
        form.send()
        msg = _('Confirmation email has been sent. '
                'Please check your inbox.')
        messages.success(request, msg)
        return redirect(settings.LOGIN_REDIRECT_URL)

    return TemplateResponse(
        request, 'registration/request_email_confirmation.jade',
        {'form': form})


def confirm_email(request, token):
    if not request.POST:
        try:
            email_confirmation_request = EmailConfirmationRequest.objects.get(
                token=token, valid_until__gte=now())
            # TODO: cronjob (celery task) to delete stale tokens
        except EmailConfirmationRequest.DoesNotExist:
            return TemplateResponse(request, 'registration/invalid_token.jade')
        user = email_confirmation_request.get_authenticated_user()
        email_confirmation_request.delete()
        auth_login(request, user)
        messages.success(request, _('You are now logged in.'))

    form = forms.SetOrRemovePasswordForm(user=request.user,
                                         data=request.POST or None)
    if form.is_valid():
        form.save()
        messages.success(request, _('Password has been successfully changed.'))
        return redirect(settings.LOGIN_REDIRECT_URL)

    return TemplateResponse(
        request, 'registration/set_password.jade', {'form': form})


def change_email(request, token):
    try:
        email_change_request = EmailChangeRequest.objects.get(
            token=token, valid_until__gte=now())
        # TODO: cronjob (celery task) to delete stale tokens
    except EmailChangeRequest.DoesNotExist:
        return TemplateResponse(request, 'registration/invalid_token.jade')

    # if another user is logged in, we need to log him out, to allow the email
    # owner confirm his identity
    if (request.user.is_authenticated() and
            request.user != email_change_request.user):
        auth_logout(request)
    if not request.user.is_authenticated():
        query = urlencode({
            'next': request.get_full_path(),
            'email': email_change_request.user.email})
        login_url = utils.url(path=settings.LOGIN_URL, query=query)
        return redirect(login_url)

    request.user.email = email_change_request.email
    request.user.save()
    email_change_request.delete()

    messages.success(request, _('Your email has been successfully changed'))
    return redirect(settings.LOGIN_REDIRECT_URL)


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
