# learn shortcuts
from django.shortcuts import get_object_or_404, render
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.views import generic
from django.utils import timezone

from django.http import HttpResponse, HttpResponseServerError
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms import AuthenticationForm
from django.shortcuts import get_object_or_404, redirect, render_to_response, RequestContext
from django.template.response import TemplateResponse
from django.contrib.messages import get_messages
from django.contrib.auth.models import User
from django.contrib import messages

from ..common.formatter import bioc_as_json, apply_bioc_annotations
from ..relationships.models import Paper, Annotation
from ..userprofile.models import UserProfile
from ..document.models import Document

from brabeion import badges
import datetime
import re

from .models import Paper, Relation, Answer
from .forms import AnswerForm
#from .tasks import parse_input
from django import forms
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


@login_required
def home(request):
    queryset_papers = Paper.objects.all()
    ctx = {
        'papers': queryset_papers,
    }
    return TemplateResponse(request, 'relationship/home.jade', ctx)

@login_required
def asdfsadfsdf(request, relationship_pk):
    model = Answer # TODO, need this here so that the questions displayed know where to POST?
    paper = get_object_or_404(Paper, pk=relationship_pk)

    chemical_anns = Annotation.objects.filter(stype='chemical').filter(paper=paper)
    # print chemical_anns, "\n", len(chemical_anns)
    disease_anns = Annotation.objects.filter(stype='disease').filter(paper=paper)
    # print disease_anns, "\n", len(disease_anns)
    relationship_pair_list = []
    relationship_pair = []
    for chemical in chemical_anns:
        for disease in disease_anns:
            # TODO there are lots of repeats here (all start stop locations are different)
            relationship_pair = [chemical, disease]
            relationship_pair_list.append(relationship_pair)

    chemical_from_relation = relationship_pair_list[0][0] # chemical from pair
    disease_from_relation = relationship_pair_list[0][1] # disease from pair
    relation = Relation.objects.first()  #TODO, make the relation here, become associated with the paper that we got.
    ctx = {'chemical': chemical_from_relation,
           'disease': disease_from_relation,
           'current_paper': paper,
           'relation': relation
           }

    # use the relationship_task.html template to work with above code given ctx
    return TemplateResponse(request, 'relationship/relationship_task.html', ctx)


#pass in relationship type pk similar to above method
def relationship_type(request, relation_id):
    # Definitely just for testing purposes. TODO fix this
    relation = Relation.objects.get(pk=relation_id)
    Answer.objects.create(relation=relation, relation_pair=relation.relation, relationship_type=request.POST['relationship_type'], user_confidence=request.POST['user_confidence'])
    return HttpResponseRedirect(reverse("relationship:home"))
