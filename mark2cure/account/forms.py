from django import forms
from django.contrib.auth.models import User
from django.forms.models import model_to_dict, fields_for_model

from mark2cure.account.models import UserProfile


class UserForm(forms.ModelForm):
    class Meta:
        model = User
        fields = ['first_name', 'last_name', 'password',
                  'email', 'username']

    def __init__(self, *args, **kwargs):
        super(UserForm, self).__init__(*args, **kwargs)
        self.fields['email'].required = True

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
        instance = kwargs.pop('instance', None)
        _fields = ('first_name', 'last_name', 'email',)
        _initial = model_to_dict(instance.user, _fields) if instance is not None else {}
        super(UserProfileForm, self).__init__(initial=_initial, instance=instance, *args, **kwargs)
        self.fields.update(fields_for_model(User, _fields))

    class Meta:
        model = UserProfile
        fields = ['timezone', 'avatar', 'email_notify',
                  'gender', 'age', 'occupation', 'education',
                  'science_education', 'country', 'referral',
                  'motivation', 'quote']
        exclude = ['user',]

    def save(self, *args, **kwargs):
        u = self.instance.user
        u.first_name = self.cleaned_data['first_name']
        u.last_name = self.cleaned_data['last_name']
        u.email = self.cleaned_data['email']
        u.save()
        profile = super(UserProfileForm, self).save(*args,**kwargs)
        return profile
