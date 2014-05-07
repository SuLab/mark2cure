from django import forms

from mark2cure.document.models import Document, Annotation, Refute, Comment


class DocumentForm(forms.ModelForm):
    class Meta:
        model = Document

        fields = ['document_id']

class AnnotationForm(forms.ModelForm):
    class Meta:
        model = Annotation

        fields = ['kind', 'text', 'start']

        # unique = ['start', 'text', 'fk_view_id']


class RefuteForm(forms.ModelForm):
    class Meta:
        model = Refute
        fields = ['message']

class CommentForm(forms.ModelForm):
    class Meta:
        model = Comment
        fields = ['message']

