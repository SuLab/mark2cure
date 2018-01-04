from django.views.decorators.http import require_http_methods
from django.shortcuts import get_object_or_404, redirect
from django.template.response import TemplateResponse
from django.http import HttpResponse

from allauth.socialaccount.models import SocialApp
from ..userprofile.models import UserProfile

from .utils.mdetect import UAgentInfo
from .forms import SupportMessageForm
from .models import Group


def get_started(request):
    '''View for directing the request into Training
    '''

    if not request.user.is_authenticated():
        return redirect('account_signup')

    return redirect('training:route')


def dashboard(request):
    '''View for serving the Dashboard.js App and issuing authenication redirects
    '''
    if not request.user.is_authenticated():
        return redirect('common:home')

    # We redirect user to the training route if their skill is not level 7
    if request.user.profile.unlocked_tasks() == 0:
        return redirect('training:route')

    uai = UAgentInfo(request.META.get('HTTP_USER_AGENT'), request.META.get('HTTP_ACCEPT'))
    return TemplateResponse(request, 'common/dashboard.html', {'mobile': uai.detectMobileLong()})


def why_mark2cure(request):
    query = UserProfile.objects.exclude(motivation='').order_by('?').values('motivation', 'user')
    return TemplateResponse(request, 'common/why-i-mark2cure.html', {'profiles': query})


def login_with_zooniverse(request):
    '''The view that allows the user to login with their Zooniverse account'''

    # The public Zooniverse app ID (used for Zoonivser login/authentication) - use the SocialApp DB record to retrieve the app ID
    zooniverseSocialApp = SocialApp.objects.get(provider='zooniverse')
    appId = zooniverseSocialApp.client_id

    ctx = {
        'zooniverse_app_id': appId,
        'zooniverse_callback_url': 'https://mark2cure.org/zooniverse-callback/'
    }
    return TemplateResponse(request, 'common/login_with_zooniverse.html', ctx)


def zooniverse_callback(request):
    '''The view that is called after the user logged-in using Zooniverse (this view receives
       the access token and sends it forward to the allauth module's callback view)'''
    ctx = {
        'zooniverse_allauth_callback_url': 'https://mark2cure.org/accounts/zooniverse/login/callback/'
    }
    return TemplateResponse(request, 'common/zooniverse_callback.html', ctx)


def ner_group_home(request, group_stub):
    '''Landing page for NER Group
    '''
    group = get_object_or_404(Group, stub=group_stub)
    return TemplateResponse(request, 'common/group_home.html', {'group': group})


@require_http_methods(['POST'])
def support(request):
    '''API for submitting SupportMessages
    '''
    form = SupportMessageForm(data=request.POST)
    if form.is_valid():
        form.save()
        return HttpResponse(200)
    return HttpResponse(500)


def home(request):
    '''Landing page for public https://mark2cure.org
    '''
    if request.user.is_authenticated():
        return redirect('common:dashboard')
    return TemplateResponse(request, 'common/landing.html')

