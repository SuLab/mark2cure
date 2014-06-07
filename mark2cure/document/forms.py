from django import forms

from mark2cure.document.models import Document, Annotation, Comment


class DocumentForm(forms.ModelForm):
    class Meta:
        model = Document
        fields = ['document_id']


class AnnotationForm(forms.ModelForm):

    '''
    def __init__(self, *args, **kwargs):
        self.view = kwargs.pop('view')
        super(AnnotationForm, self).__init__(*args, **kwargs)
        self.fields['view_id'].initial = self.view.pk
    '''

    class Meta:
        model = Annotation
        fields = ['kind', 'text', 'start']


class CommentForm(forms.ModelForm):
    class Meta:
        model = Comment
        fields = ['message']

