from django import forms
from django.contrib.auth.models import User

from mark2cure.common.models import Message


class MessageForm(forms.ModelForm):
    class Meta:
        model = Message
        fields = ['message']
