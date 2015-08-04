import datetime
import re

# learn shortcuts
from django.shortcuts import get_object_or_404, render
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.views import generic
from django.utils import timezone


from .models import Paper
#from .tasks import parse_input

import os
# TODO check imports
### add here if not completed TODO (remove from list if user already answered)
"""
Toby's parse_input creates new Paper objects using all of the
inputs to the Paper class.
"""
# TODO, views for this handled in common?  


class DetailView(generic.DetailView):
    model = Paper
    template_name = 'relationships/detail.html'

    def get_queryset(self):
        return Paper.objects.all()

class ResultsView(generic.DetailView):
    model = Paper
    template_name = 'relationships/results.html'
