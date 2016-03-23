from django import forms
from .models import EntityRecognitionAnnotation


class EntityRecognitionAnnotationForm(forms.ModelForm):
    class Meta:
        model = EntityRecognitionAnnotation
        fields = ['text', 'start', 'type']

    def clean_type(self):
        data = self.cleaned_data['type']
        if data.isdigit():
            data = EntityRecognitionAnnotation.ANNOTATION_TYPE_CHOICE[int(data)]
        return data
