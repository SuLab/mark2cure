from django.contrib import messages
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse, reverse_lazy
from django.contrib.auth.decorators import login_required
from django.views.generic.base import TemplateView
from django.views.generic.edit import FormView

from allauth.account import app_settings as account_settings
from allauth.account.views import (AjaxCapableProcessFormViewMixin,
                             CloseableSignupMixin,
                             RedirectAuthenticatedUserMixin)
from allauth.account.adapter import get_adapter as get_account_adapter
from allauth.utils import get_form_class, get_current_site

from allauth.socialaccount.adapter import get_adapter
from allauth.socialaccount.models import SocialLogin
from allauth.socialaccount.forms import DisconnectForm, SignupForm
from allauth.socialaccount import helpers
from allauth.socialaccount import app_settings


class SignupView(RedirectAuthenticatedUserMixin, CloseableSignupMixin,
                 AjaxCapableProcessFormViewMixin, FormView):
    form_class = SignupForm
    template_name = (
        'socialaccount/signup.' + account_settings.TEMPLATE_EXTENSION)

    def get_form_class(self):
        return get_form_class(app_settings.FORMS,
                              'signup',
                              self.form_class)

    def dispatch(self, request, *args, **kwargs):
        if not request.session.get('initial_training', False):
            return redirect('common:home')

        self.sociallogin = None
        data = request.session.get('socialaccount_sociallogin')
        if data:
            self.sociallogin = SocialLogin.deserialize(data)
        if not self.sociallogin:
            return HttpResponseRedirect(reverse('account_login'))
        return super(SignupView, self).dispatch(request, *args, **kwargs)

    def is_open(self):
        return get_adapter(self.request).is_open_for_signup(
            self.request,
            self.sociallogin)

    def get_form_kwargs(self):
        ret = super(SignupView, self).get_form_kwargs()
        ret['sociallogin'] = self.sociallogin
        return ret

    def form_valid(self, form):
        form.save(self.request)

        # After loggin them in, assign the first Level training so we know where to route them
        Level.objects.create(user=self.user, task_type=self.request.session.get('initial_training'), level=3, created=timezone.now())

        return helpers.complete_social_signup(self.request,
                                              self.sociallogin)

    def get_context_data(self, **kwargs):
        ret = super(SignupView, self).get_context_data(**kwargs)
        ret.update(dict(site=get_current_site(self.request),
                        account=self.sociallogin.account))
        return ret

    def get_authenticated_redirect_url(self):

        return reverse(connections)

signup = SignupView.as_view()

