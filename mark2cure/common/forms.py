from django import forms
from .models import SupportMessage


class SupportMessageForm(forms.ModelForm):
    class Meta:
        model = SupportMessage
        fields = ['user', 'text', 'referral']

