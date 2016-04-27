from django.core.urlresolvers import reverse
from django.views.generic.edit import FormView
from django.shortcuts import redirect
from django.views.decorators.debug import sensitive_post_parameters
from django.utils.decorators import method_decorator

from allauth.exceptions import ImmediateHttpResponse
from allauth.account.utils import (get_next_redirect_url, complete_signup,
                                        get_login_redirect_url, perform_login,
                                        passthrough_next_redirect_url, url_str_to_user_pk,
                                        logout_on_password_change)

from allauth.utils import get_form_class, get_request_param


from allauth.account.forms import (
    AddEmailForm, ChangePasswordForm,
    LoginForm, ResetPasswordKeyForm,
    ResetPasswordForm, SetPasswordForm, SignupForm, UserTokenForm)

from allauth.account import app_settings
from allauth.account.adapter import get_adapter

from ...task.models import Level
from django.utils import timezone


sensitive_post_parameters_m = method_decorator(sensitive_post_parameters('password', 'password1', 'password2'))

def _ajax_response(request, response, form=None):
    if request.is_ajax():
        if (isinstance(response, HttpResponseRedirect) or isinstance(response, HttpResponsePermanentRedirect)):
            redirect_to = response['Location']
        else:
            redirect_to = None
        response = get_adapter(request).ajax_response(request, response, form=form, redirect_to=redirect_to)
    return response


class RedirectAuthenticatedUserMixin(object):
    def dispatch(self, request, *args, **kwargs):
        # WORKAROUND: https://code.djangoproject.com/ticket/19316
        self.request = request
        # (end WORKAROUND)
        if request.user.is_authenticated() and app_settings.AUTHENTICATED_LOGIN_REDIRECTS:
            redirect_to = self.get_authenticated_redirect_url()
            response = HttpResponseRedirect(redirect_to)
            return _ajax_response(request, response)
        else:
            response = super(RedirectAuthenticatedUserMixin, self).dispatch(request, *args, **kwargs)
        return response

    def get_authenticated_redirect_url(self):
        redirect_field_name = self.redirect_field_name
        return get_login_redirect_url(self.request, url=self.get_success_url(), redirect_field_name=redirect_field_name)


class AjaxCapableProcessFormViewMixin(object):

    def post(self, request, *args, **kwargs):
        form_class = self.get_form_class()
        form = self.get_form(form_class)
        if form.is_valid():
            response = self.form_valid(form)
        else:
            response = self.form_invalid(form)
        return _ajax_response(self.request, response, form=form)


class SignupView(RedirectAuthenticatedUserMixin, AjaxCapableProcessFormViewMixin, FormView):
    template_name = "account/signup." + app_settings.TEMPLATE_EXTENSION
    form_class = SignupForm
    redirect_field_name = "next"
    success_url = None

    @sensitive_post_parameters_m
    def dispatch(self, request, *args, **kwargs):
        if not request.session.get('initial_training', False):
            return redirect('common:home')

        return super(SignupView, self).dispatch(request, *args, **kwargs)

    def get_form_class(self):
        return get_form_class(app_settings.FORMS, 'signup', self.form_class)

    def get_success_url(self):
        # After loggin them in, assign the first Level training so we know where to route them
        Level.objects.create(user=self.user, task_type=self.request.session.get('initial_training'), level=3, created=timezone.now())

        # Explicitly passed ?next= URL takes precedence
        ret = (
            get_next_redirect_url(
            self.request,
            self.redirect_field_name) or self.success_url)
        return ret

    def form_valid(self, form):
        # By assigning the User to a property on the view, we allow subclasses
        # of SignupView to access the newly created User instance
        self.user = form.save(self.request)
        return complete_signup( self.request, self.user,
                                app_settings.EMAIL_VERIFICATION,
                                self.get_success_url())

    # def post(self, request, *args, **kwargs):
    #    print '> signup post empty', self

    def get_context_data(self, **kwargs):
        ret = super(SignupView, self).get_context_data(**kwargs)
        form = ret['form']
        form.fields["email"].initial = self.request.session.get('account_verified_email')
        login_url = passthrough_next_redirect_url(self.request, reverse("account_login"), self.redirect_field_name)

        redirect_field_name = self.redirect_field_name
        redirect_field_value = get_request_param(self.request, redirect_field_name)

        ret.update({"login_url": login_url, "redirect_field_name": redirect_field_name, "redirect_field_value": redirect_field_value})
        return ret

signup = SignupView.as_view()

