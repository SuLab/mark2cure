from django import forms
from .models import EntityRecognitionAnnotation


class EntityRecognitionAnnotationForm(forms.ModelForm):
    '''
        Responsible for saving POST'd data from the YPet
        library.

        - The start position that is given is relative to the
        passage it's contained within, NOT the entire document
    '''

    class Meta:
        model = EntityRecognitionAnnotation
        fields = ['text', 'start', 'type']

    def clean_type(self):
        data = self.cleaned_data['type']
        if data.isdigit():
            data = EntityRecognitionAnnotation.ANNOTATION_TYPE_CHOICE[int(data)]
        return data
