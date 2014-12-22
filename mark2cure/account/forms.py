from django import forms
from django.contrib.auth.models import User

from mark2cure.account.models import UserProfile
from django.contrib.auth.forms import UserCreationForm


class UserCreateForm(UserCreationForm):
    email = forms.EmailField(required=True)

    class Meta:
        model = User
        fields = ("username", "email", "password1", "password2")

    def save(self, commit=True):
        user = super(UserCreateForm, self).save(commit=False)
        user.email = self.cleaned_data["email"]
        if commit:
            user.save()
            return user


class UserNameChangeForm(forms.ModelForm):
    class Meta:
        model = User
        fields = ('first_name', 'last_name')


class UserProfileForm(forms.ModelForm):
    class Meta:
        model = UserProfile
        fields = ['timezone', 'avatar', 'email_notify',
                  'gender', 'age', 'occupation', 'education',
                  'science_education', 'country', 'referral',
                  'motivation', 'quote']
        exclude = ['user', ]

