from __future__ import unicode_literals

from django import forms
from django.contrib.auth.models import User
from django.utils.translation import ugettext

from django.core.mail import send_mail

from django.contrib.auth.forms import UserCreationForm, SetPasswordForm

from .models import EmailConfirmationRequest, EmailChangeRequest


class UserCreateForm(UserCreationForm):
    email = forms.EmailField(required=True)

    class Meta:
        model = User
        fields = ('username', 'email', 'password1', 'password2')

    def save(self, commit=True):
        user = super(UserCreateForm, self).save(commit=False)
        user.email = self.cleaned_data['email']
        if commit:
            user.save()
            return user


class UserNameChangeForm(forms.ModelForm):
    class Meta:
        model = User
        fields = ('first_name', 'last_name')


class SetOrRemovePasswordForm(SetPasswordForm):

    def __init__(self, *args, **kwargs):
        super(SetOrRemovePasswordForm, self).__init__(*args, **kwargs)
        if 'new_password1' not in self.data.keys():
            self.fields['new_password1'].required = False
            self.fields['new_password2'].required = False

    def save(self, commit=True):
        if self.cleaned_data.get('new_password1'):
            return super(SetOrRemovePasswordForm, self).save(commit)
        else:
            self.user.set_unusable_password()
        return self.user


class RequestEmailConfirmationForm(forms.Form):

    email = forms.EmailField()

    template = 'registration/emails/confirm_email.txt'

    def __init__(self, local_host=None, data=None):
        self.local_host = local_host
        super(RequestEmailConfirmationForm, self).__init__(data)

    def send(self):
        email = self.cleaned_data['email']
        request = self.create_request_instance()
        confirmation_url = self.local_host + request.get_confirmation_url()
        context = {'confirmation_url': confirmation_url}
        send_mail(email, self.template, context)

    def create_request_instance(self):
        email = self.cleaned_data['email']
        EmailConfirmationRequest.objects.filter(email=email).delete()
        return EmailConfirmationRequest.objects.create(
            email=self.cleaned_data['email'])


class RequestEmailChangeForm(RequestEmailConfirmationForm):

    template = 'registration/emails/change_email.txt'

    def __init__(self, user=None, *args, **kwargs):
        self.user = user
        super(RequestEmailChangeForm, self).__init__(*args, **kwargs)

    def clean_email(self):
        email = self.cleaned_data['email']
        if User.objects.filter(email=email).exists():
            raise forms.ValidationError(
                ugettext('Account with this email already exists'))
        return self.cleaned_data['email']

    def create_request_instance(self):
        EmailChangeRequest.objects.filter(user=self.user).delete()
        return EmailChangeRequest.objects.create(
            email=self.cleaned_data['email'], user=self.user)

