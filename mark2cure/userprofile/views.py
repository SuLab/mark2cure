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

from .forms import UserNameChangeForm


@login_required
def settings(request):
    ctx = {'menu': 4}
    return TemplateResponse(request, 'userprofile/myinfo.jade', ctx)


@login_required
def settings_edit(request):
    user_name_form = UserNameChangeForm(instance=request.user, data=request.POST or None)

    if request.POST and user_name_form.is_valid():
        user_name_form.save()
        return redirect('userprofile:settings')

    ctx = { 'user_name_form': user_name_form,
            'menu': 4}
    return TemplateResponse(request, 'userprofile/myinfo-edit.jade', ctx)


@login_required
def paypal(request):
    ctx = {'menu': 4}
    return TemplateResponse(request, 'userprofile/paypalinfo.jade', ctx)


@login_required
def receipts(request):
    ctx = {'menu': 4}
    return TemplateResponse(request, 'userprofile/myinfo.jade', ctx)

