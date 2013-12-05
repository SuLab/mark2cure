'''
django model forms allow you to
easily tie models to views

e.g. for doc we should have something like
'''

from django import forms
from django.core.exceptions import ObjectDoesNotExist

from mark2cure.document.models import Document, Annotation, View
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

    # def __init__(self, *args, **kwargs):
    #     v = kwargs.pop('userview', None)
    #     print "\n\n View is here: ", v

    #     super(AnnotationForm, self).__init__(*args, **kwargs)
    #     # forms.ModelForm.__init__(self, *args, **kwargs)
    #     self.view = v


