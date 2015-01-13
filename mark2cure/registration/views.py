try:
    from urllib.parse import urlencode
except ImportError:
    from urllib import urlencode

from django.conf import settings
from django.contrib import messages
from django.contrib.auth import authenticate, login as auth_login, logout as auth_logout

from django.contrib.auth.decorators import login_required

from django.contrib.auth.views import (
    login as django_login_view, password_change)
from django.core.urlresolvers import reverse
from django.shortcuts import get_object_or_404, render_to_response, redirect
from django.template.response import TemplateResponse
from django.utils import timezone
from django.utils.encoding import smart_text
from django.utils.translation import ugettext_lazy as _

from django.views.decorators.http import require_http_methods
from django.views.decorators.csrf import csrf_exempt
from django.contrib.auth import get_user_model

from django.template import RequestContext

from . import forms
from .models import EmailConfirmationRequest, EmailChangeRequest
from .forms import UserCreateForm

from ecleanit.userprofile.forms import UserProfileForm

from ecleanit.cleanit.models import Property, Cleaner, Job
from ecleanit.cleanit.forms import PropertyForm, PropertySchedulePerenceForm, CleanerForm, JobForm
from ecleanit.userprofile.models import UserProfile
from ecleanit.userprofile.forms import UserProfileForm

from . import utils

now = timezone.now


def user_creation(request):
    if request.user.is_authenticated():
        return redirect('registration:setup-property')

    user_create_form = UserCreateForm(data=request.POST or None)

    if request.POST and user_create_form.is_valid():
        user = user_create_form.save()
        user = authenticate(
                username=request.POST['username'],
                password=request.POST['password1'])
        auth_login(request, user)

        return redirect('registration:setup-property')


    return TemplateResponse(request, 'registration/register-step-1.jade',
                            {'form': user_create_form})


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


def change_password(request):
    return password_change(
        request, template_name='registration/change_password.jade',
        post_change_redirect=reverse('userprofile:settings'))
