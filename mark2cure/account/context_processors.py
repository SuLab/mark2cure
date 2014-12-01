# -*- coding: utf-8 -*-
from mark2cure.account.forms import UserForm, UserProfileForm
from django.contrib.auth.forms import AuthenticationForm

def inject_signup_forms(request):
    if not request.user.is_authenticated():
        return {'create_form': UserForm(),
                'login_form': AuthenticationForm()}
    else:
        return {}
