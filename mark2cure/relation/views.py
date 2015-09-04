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
import json
import random

from .models import Paper, Relation, Answer
from django import forms
import os


# TODO (remove this, if I am going to use the JQUERY vertical menu)
my_json ={'associations': {
            "chemical_disease": {
                "exacerbates": [],
                "increases risk of": [],
                "may cause": [],
                "prevents": [],
                "treats": {
                    'treats symptoms of': [],
                    'treats cause of': []
                    }
                },
            "chemical_gene": {
                'binds to': [],
                'changes location of': [],
                'metabolizes': [],
                'transports': [],
                'not sure about relation': []
                },
            'gene_disease': {
                'altered expression': [],
                'altered regulation leads to': [],
                'mutation association': [],
                'post translational modification association':[],
                'not sure about relation': []
                }

            }
}

json_string = json.dumps(my_json)
parsed_json = json.loads(json_string)


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

    # TODO, this might be a "weird" chemical name so take the majority name?
    chemical_from_relation = Annotation.objects.filter(uid=relation.chemical_id).filter(paper=current_paper.id)[0].text
    disease_from_relation = Annotation.objects.filter(uid=relation.disease_id).filter(paper=current_paper.id)[0].text

    formatted_abstract = re.sub(chemical_from_relation, '<font color="#E65CE6"><b>'+chemical_from_relation+'</b></font>', current_paper.abstract, flags=re.I)
    formatted_abstract = re.sub(disease_from_relation, '<font color="#0099FF"><b>'+disease_from_relation+'</b></font>', formatted_abstract, flags=re.I)

    chemical_from_relation = '<font color="#E65CE6"><b>'+ chemical_from_relation + '</b></font>'
    disease_from_relation = '<font color="#0099FF"><b>'+ disease_from_relation+ '</b></font>'

    # TODO: want something similar to this:
    # paper = relation_task.papers().first()
    # sentences = relation_task.get_sentences()

    question_list1 = parsed_json['associations']['chemical_disease']
    json_object = parsed_json['associations']['chemical_disease']['treats']
    print json_object

    level = 2  #(TODO) remove this "level": here as the new button name, this will be removed by widget.

    # print parsed_json['associations']['chemical_disease']['treats']['treats cause']

    ctx = {'chemical' : chemical_from_relation,
           'disease' : disease_from_relation,
           'json_object': json_object,
           'level': level,
           'question_list': question_list1,
           'current_paper': formatted_abstract,
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
