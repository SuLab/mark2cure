from django import forms

from mark2cure.document.models import Document, Annotation


class DocumentForm(forms.ModelForm):
    class Meta:
        model = Document
        fields = ['document_id']


class AnnotationForm(forms.ModelForm):
    class Meta:
        model = Annotation
        fields = ['text', 'start', 'type']

    def clean_type(self):
        data = self.cleaned_data['type']
        if data.isdigit():
            data = Annotation.ANNOTATION_TYPE_CHOICE[int(data)]
        return data
