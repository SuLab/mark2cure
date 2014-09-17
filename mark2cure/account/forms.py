from django import forms
from django.contrib.auth.models import User

from mark2cure.account.models import UserProfile


class UserForm(forms.ModelForm):

    class Meta:
        model = User
        fields = ['first_name', 'last_name', 'password',
                  'email', 'username']


    def save(self): # create new user
        new_user = User.objects.create_user(
            username=self.cleaned_data['username'],
            password=self.cleaned_data['password'],
            email=self.cleaned_data['email'])

        new_user.first_name = self.cleaned_data['first_name']
        new_user.last_name = self.cleaned_data['last_name']
        new_user.save()

        return new_user


class UserProfileForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        super(UserProfileForm, self).__init__(*args, **kwargs)

        # for key in self.fields:
        #     self.fields[key].required = True

    class Meta:
        model = UserProfile
        fields = ['timezone', 'avatar', 'email_notify',
                  'gender', 'age', 'occupation', 'education',
                  'science_education', 'country', 'referral',
                  'motivation', 'quote']
