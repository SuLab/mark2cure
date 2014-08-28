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

from mark2cure.account.models import UserProfile, Message
from mark2cure.account.forms import UserForm, ProfileForm

import datetime


def reset_thanks(request):
    return render_to_response('account/reset-thanks.jade', {}, context_instance=RequestContext(request))


@login_required
def settings(request):
    user = request.user
    profileForm = ProfileForm(instance=user.profile)
    passwordChangeForm = PasswordChangeForm(user)

    print passwordChangeForm

    return render_to_response('account/settings.jade',
            {'passwordChangeForm': passwordChangeForm,
             'profileForm': profileForm},
              context_instance=RequestContext(request))


@login_required
@require_http_methods(['POST'])
def update_profile(request, profile_id):
    profile = get_object_or_404(UserProfile, pk=profile_id)
    form = ProfileForm(request.POST, instance=profile)

    if form.is_valid():
        # print form
        profile = form.save()

    return redirect('/account/')


def create(request):
    if request.method == 'POST':
        form = UserForm(request.POST)
        if form.is_valid():
            user = form.save()
            profile = user.profile
            #profile.created_by = request.user
            profile.save()

            user = authenticate(username=request.POST['username'], password=request.POST['password'])
            login(request, user)
            return redirect('mark2cure.common.views.quest_read')

    else:
        form = UserForm()

    return render_to_response('account/create.jade', {'form': form}, context_instance=RequestContext(request))


@login_required
@require_http_methods(['POST'])
def inactivate(request, user_id):
    if not (request.user.is_staff or request.user.is_superuser):
        return HttpResponse('Unauthorized', status=401)

    user = get_object_or_404(User, pk=user_id)

    if request.user.is_superuser or user.profile.created_by.pk == request.user.pk:
        user.is_active = False
        user.save()
        return redirect('/account/')

    else:
        return HttpResponse('Unauthorized', status=401)


@login_required
@require_http_methods(['POST'])
def activate(request, user_id):
    if not (request.user.is_staff or request.user.is_superuser):
        return HttpResponse('Unauthorized', status=401)

    user = get_object_or_404(User, pk=user_id)

    if request.user.is_superuser or user.profile.created_by.pk == request.user.pk:
        user.is_active = True
        user.save()
        return redirect('/account/')

    else:
        return HttpResponse('Unauthorized', status=401)


@login_required
@require_http_methods(['POST'])
def delete(request, user_id):
    if not (request.user.is_staff or request.user.is_superuser):
        return HttpResponse('Unauthorized', status=401)

    user = get_object_or_404(User, pk=user_id)

    if request.user.is_superuser or user.profile.created_by.pk == request.user.pk:
        user.delete()
        return redirect('/account/')

    else:
        return HttpResponse('Unauthorized', status=401)



@login_required
@require_http_methods(['POST'])
def staffify(request, user_id):
    if not (request.user.is_staff or request.user.is_superuser):
        return HttpResponse('Unauthorized', status=401)

    user = get_object_or_404(User, pk=user_id)

    if request.user.is_superuser or user.profile.created_by.pk == request.user.pk:
        user.is_staff = True
        user.save()
        return redirect('/account/')

    else:
        return HttpResponse('Unauthorized', status=401)


@login_required
@require_http_methods(['POST'])
def destaffify(request, user_id):
    if not (request.user.is_staff or request.user.is_superuser):
        return HttpResponse('Unauthorized', status=401)

    user = get_object_or_404(User, pk=user_id)

    if request.user.is_superuser or user.profile.created_by.pk == request.user.pk:
        user.is_staff = False
        user.save()
        return redirect('/account/')

    else:
        return HttpResponse('Unauthorized', status=401)

