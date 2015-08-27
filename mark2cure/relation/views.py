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
#from .forms import AnswerForm
#from .tasks import parse_input
from django import forms
import os


@login_required
def home(request):
    queryset_papers = Paper.objects.all()

    exclude_list = []
    for i in queryset_papers:
        user_paper_relations = Relation.objects.filter(paper=i.id)
        flag = 0
        for j in user_paper_relations:
            user_paper_answers = Answer.objects.all().filter(relation_pair=j.relation).filter(username=request.user)
            if len(user_paper_answers) == 1:
                flag += 1
            if flag == len(user_paper_relations):
                exclude_list.append(i)

    exclude_these = [i.id for i in exclude_list]
    queryset_papers = Paper.objects.all().exclude(id__in=exclude_these)

    ctx = {
        'papers': queryset_papers,
    }
    return TemplateResponse(request, 'relation/home.jade', ctx)


@login_required
def relation(request, paper_pk):
    model = Answer # TODO, need this here so that the questions displayed know where to POST?
    current_paper = get_object_or_404(Paper, pk=paper_pk)
    relations = Relation.objects.filter(paper=current_paper.id)


    for i in relations:
        relation = i
        # relations that user has already completed
        relation_specific_answers = Answer.objects.filter(username=request.user).filter(relation_pair=relation.relation)
        if not relation_specific_answers:
            break

    # TODO, this might be a "weird" chemical so take the tie breakers
    chemical_from_relation = Annotation.objects.filter(uid=relation.chemical_id).filter(paper=current_paper.id)[0].text
    disease_from_relation = Annotation.objects.filter(uid=relation.disease_id).filter(paper=current_paper.id)[0].text

    # TODO: want something similar to this:
    # paper = relation_task.papers().first()
    # sentences = relation_task.get_sentences()

    ctx = {'chemical': chemical_from_relation,
           'disease': disease_from_relation,
           'current_paper': current_paper,
           'relation': relation
           }

    # use the relation.html template to work with above code given ctx
    return TemplateResponse(request, 'relation/relation.jade', ctx)


#pass in relation similar to above method
def results(request, relation_id):
    # Definitely just for testing purposes. TODO fix this
    relation = Relation.objects.get(pk=relation_id)
    current_answer = Answer.objects.create(relation=relation, relation_pair=relation.relation, relation_type=request.POST['relation_type'], user_confidence='C1', username=request.user)

    #TODO dictionary for max:
    user_work = {'relation': relation,
                 'relation_pair': relation.relation,  # not really needed but nice to see for now in th DB
                 'relation_type': request.POST['relation_type'],
                 'username': request.user
    }
    print user_work
    #this is how to display the field and not the short choice abbreviation

    ctx = {'relation_type': current_answer.get_relation_type_display()
    }
    return TemplateResponse(request, 'relation/results.jade', ctx)
    #return HttpResponseRedirect(reverse("relation:home"))
