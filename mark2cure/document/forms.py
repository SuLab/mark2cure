from django import forms
from django.core.exceptions import ObjectDoesNotExist

from mark2cure.document.models import Document, Annotation, View, Refute, Comment
from mark2cure.common.utils import get_timezone_offset

import re


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

