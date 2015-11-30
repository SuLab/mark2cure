from django import forms

from mark2cure.document.models import Annotation


class AnnotationForm(forms.ModelForm):
    class Meta:
        model = Annotation
        fields = ['text', 'start', 'type']

    def clean_type(self):
        data = self.cleaned_data['type']
        if data.isdigit():
            data = Annotation.ANNOTATION_TYPE_CHOICE[int(data)]
        return data
