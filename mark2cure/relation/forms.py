from django import forms

from .models import Answer

class AnswerForm(forms.Form):
    class Meta:
        model = Answer
        """
        fields = ('relation', 'relation_type', 'user_confidence')
        relation = forms.TextField(blank=False)
        relation_type = forms.TextField(default='relation_type')
        user_confidence = forms.TextField(default='user_confidence')
        #relation_type = forms.CharField(label='answer_1', max_length=100)
        """
