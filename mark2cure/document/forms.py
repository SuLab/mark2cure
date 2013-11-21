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

