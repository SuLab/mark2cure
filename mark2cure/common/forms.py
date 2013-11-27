from django import forms
from django.contrib.auth.models import User

class UserForm(forms.ModelForm):

    class Meta:
        model = User
        fields = ['email']


    def save(self): # signup a new user
        new_user = User.objects.create_user(username = self.cleaned_data['email'],
                                            email = self.cleaned_data['email'])
        new_user.save()
        return new_user
