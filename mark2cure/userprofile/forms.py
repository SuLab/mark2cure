from django import forms
from django.utils.translation import ugettext as _

from django.contrib.auth.models import User
#from .models import Address, AddressBook,
from .models import UserProfile


class UserProfileForm(forms.ModelForm):
    class Meta:
        model = UserProfile
        fields = ['timezone', 'avatar', 'email_notify',
                  'gender', 'age', 'occupation', 'education',
                  'science_education', 'country', 'referral',
                  'motivation', 'quote']
        exclude = ['user', ]

