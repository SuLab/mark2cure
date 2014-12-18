# -*- coding: utf-8 -*-
from django.contrib.auth.forms import AuthenticationForm

def inject_signup_forms(request):
    if not request.user.is_authenticated():
        return {'login_form': AuthenticationForm()}
    else:
        return {}
