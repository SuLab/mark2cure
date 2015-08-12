from django import forms

from .models import Answer

class AnswerForm(forms.Form):
    class Meta:
        model = Answer
        """
        fields = ('relation', 'relationship_type', 'user_confidence')
        relation = forms.TextField(blank=False)
        relationship_type = forms.TextField(default='relationship_type')
        user_confidence = forms.TextField(default='user_confidence')
        #relationship_type = forms.CharField(label='answer_1', max_length=100)
        """
