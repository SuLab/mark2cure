from django import forms
from django.contrib.auth.models import User
from mark2cure.account.models import UserProfile


from mark2cure.common.models import Message

class ProfileSurveyForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        super(ProfileSurveyForm, self).__init__(*args, **kwargs)

        for key in self.fields:
            self.fields[key].required = True

    class Meta:
        model = UserProfile
        fields = ['country', 'gender', 'age', 'occupation', 'education', 'science_education', 'motivation']


class MessageForm(forms.ModelForm):
    class Meta:
        model = Message
        fields = ['message']
