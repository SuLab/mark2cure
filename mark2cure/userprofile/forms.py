from django import forms
from .models import UserProfile, Team


class UserProfileForm(forms.ModelForm):
    class Meta:
        model = UserProfile
        fields = ['timezone', 'avatar', 'team', 'email_notify',
                  'gender', 'age', 'occupation', 'education',
                  'science_education', 'country', 'referral',
                  'motivation', 'quote']
        exclude = ['user', ]
        widgets = {
            'motivation': forms.Textarea(attrs={'rows': 2}),
            'referral': forms.Textarea(attrs={'rows': 2}),
            'quote': forms.Textarea(attrs={'rows': 4}),
        }


class TeamForm(forms.ModelForm):
    class Meta:
        model = Team
        fields = ['name', ]
