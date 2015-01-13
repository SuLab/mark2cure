from django import forms
from django.utils.translation import ugettext as _

from django.contrib.auth.models import User
#from .models import Address, AddressBook,
from .models import UserProfile


class UserProfileForm(forms.ModelForm):

    class Meta:
        model = UserProfile
        exclude = ['user']


class UserNameChangeForm(forms.ModelForm):

    class Meta:
        model = User
        fields = ('first_name', 'last_name')

