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
from ..relation.models import Paper, Annotation
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


@login_required
def home(request):
    queryset_papers = Paper.objects.all()
    ctx = {
        'papers': queryset_papers,
    }
    return TemplateResponse(request, 'relation/home.jade', ctx)


@login_required
def relation(request, paper_pk):
    model = Answer # TODO, need this here so that the questions displayed know where to POST?
    paper = get_object_or_404(Paper, pk=paper_pk)

    print request.user

    # paper = relation_task.papers().first()
    # sentences = relatinshiop_task.get_sentences()

    chemical_anns = Annotation.objects.filter(stype='chemical').filter(paper=paper)
    # print chemical_anns, "\n", len(chemical_anns)
    disease_anns = Annotation.objects.filter(stype='disease').filter(paper=paper)
    # print disease_anns, "\n", len(disease_anns)
    relation_pair_list = []
    relation_pair = []
    for chemical in chemical_anns:
        for disease in disease_anns:
            # TODO there are lots of repeats here (all start stop locations are different)
            relation_pair = [chemical, disease]
            relation_pair_list.append(relation_pair)

    chemical_from_relation = relation_pair_list[0][0] # chemical from pair
    disease_from_relation = relation_pair_list[0][1] # disease from pair
    relation = Relation.objects.first()  #TODO, make the relation here, become associated with the paper that we got.
    ctx = {'chemical': chemical_from_relation,
           'disease': disease_from_relation,
           'current_paper': paper,
           'relation': relation
           }

    # use the relation.html template to work with above code given ctx
    # TODO make these jade, not HTML.
    return TemplateResponse(request, 'relation/relation.html', ctx)


#pass in relation type pk similar to above method
def results(request, relation_id):
    # Definitely just for testing purposes. TODO fix this
    relation = Relation.objects.get(pk=relation_id)
    Answer.objects.create(relation=relation, relation_pair=relation.relation, relation_type=request.POST['relation_type'], user_confidence=request.POST['user_confidence'])

    ctx = {}
    return TemplateResponse(request, 'relation/results.jade', ctx)
    #return HttpResponseRedirect(reverse("relation:home"))
