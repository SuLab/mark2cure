from django import forms
from .models import EntityRecognitionAnnotation


class EntityRecognitionAnnotationForm(forms.ModelForm):
    '''
        Responsible for saving POST'd data from the YPet
        library.

        - The start position that is given is relative to the
        passage it's contained within, NOT the entire document
    '''
    type_id = forms.IntegerField()

    class Meta:
        model = EntityRecognitionAnnotation
        fields = ['text', 'start', 'type']

    def clean_type(self):
        # (TODO) Cast or modified this posted *_id attribute to the model
        # data = self.cleaned_data['type_id']
        data = self.data.get('type_id')
        if data.isdigit():
            data = EntityRecognitionAnnotation.ANNOTATION_TYPE_CHOICE[int(data)]
        return data
