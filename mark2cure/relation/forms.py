from django import forms
from .models import Answer

class AnswerForm(forms.ModelForm):
    class Meta:
        model = Answer
        fields = ['relation', 'relation_pair', 'relation_type', 'username']
