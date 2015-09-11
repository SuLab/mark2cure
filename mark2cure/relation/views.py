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
from ..common.bioc import BioCReader
from ..relation.models import Paper, Annotation
from ..userprofile.models import UserProfile
from ..document.models import Document, Section

from brabeion import badges
import datetime
import re
import json
import random

from .models import Answer, Relation
from ..document.models import Document, Pubtator
from .forms import AnswerForm
# from django import forms
import os


# TODO (remove this, if I am going to use the JQUERY vertical menu)


def make_annotation_lists_from_current_document(current_document):
    """
    TODO: should be removed from views and possibly added to models or tasks?
    Input is one document. Output is three concept lists in order: disease, gene,
    chemical. Also output is the pubtator_anns.
    """
    concept0_list = []
    concept1_list = []
    concept2_list = []
    pubtator_anns = ''
    if current_document.valid_pubtator():
        pubtator_anns = current_document.as_bioc_with_pubtator_annotations()
        num_passages = len(pubtator_anns.passages)
        for i in range(0, num_passages):
            len_annotations = len(pubtator_anns.passages[i].annotations)
            for k in range(0, len_annotations):
                concept_text = pubtator_anns.passages[i].annotations[k].text
                concept_type = pubtator_anns.passages[i].annotations[k].infons['type']
                if concept_type == "0":
                    concept0_list.append(concept_text)
                if concept_type == "1":
                    concept1_list.append(concept_text)
                if concept_type == "2":
                    concept2_list.append(concept_text)
    return concept0_list, concept1_list, concept2_list, pubtator_anns



@login_required
def home(request):
    # how many documents to display on the home page
    queryset_papers = Document.objects.all()[0:10]

    filtered_queryset_papers = []
    for current_document in queryset_papers:
        current_document = Document.objects.get(pk=current_document.pk)
        concept0_list, concept1_list, concept2_list, pubtator_anns = make_annotation_lists_from_current_document(current_document)
        if len(concept0_list) != 0 and len(concept2_list) != 0:
            filtered_queryset_papers.append(current_document)

    queryset_papers = filtered_queryset_papers


    # for i in queryset_papers:
    #     # print i.id
    #     # print i.document_id
    #     my_pub = Pubtator.objects.filter(document=i)[0]



    # exclude if user already answered
    # exclude_list = []
    # for i in queryset_papers:
    #     user_paper_relations = Relation.objects.filter(paper=i.id)
    #     print user_paper_relations
    #     flag = 0
    #     for j in user_paper_relations:
    #         user_paper_answers = Answer.objects.all().filter(relation_pair=j.relation).filter(username=request.user)
    #         if len(user_paper_answers) == 1:
    #             flag += 1
    #         if flag == len(user_paper_relations):
    #             exclude_list.append(i)
    #
    # exclude_these = [i.id for i in exclude_list]
    # # queryset_papers = Paper.objects.all().exclude(id__in=exclude_these)

    ctx = {
        'documents': queryset_papers,
    }
    return TemplateResponse(request, 'relation/home.jade', ctx)


@login_required
def relation(request, document_pk):
    print "\n\n\n\n\n\n\n"
    model = Answer # TODO, need this here so that the questions displayed know where to POST?
    current_document = get_object_or_404(Document, pk=document_pk)
    # relations = Relation.objects.filter(paper=current_document.id)
    # print title, "title\n\n"
    # print abstract, 'abstract\n\n'

    concept0_list, concept1_list, concept2_list, pubtator_anns = make_annotation_lists_from_current_document(current_document)
    # disease, gene, chemical
    print concept0_list
    print concept1_list
    print concept2_list


    #
    # for relation in relations:
    #     # relations that user has already completed
    #     relation_specific_answers = Answer.objects.filter(username=request.user).filter(relation_pair=relation.relation)
    #     if not relation_specific_answers:
    #         break

    # # TODO, this might be a "weird" chemical name so take the majority name?
    # chemical_from_relation = Annotation.objects.filter(uid=relation.chemical_id).filter(paper=current_paper.id)[0].text
    # disease_from_relation = Annotation.objects.filter(uid=relation.disease_id).filter(paper=current_paper.id)[0].text

    chemical = concept2_list[0]
    disease = concept0_list[0]

    formatted_abstract = re.sub(r'\b'+chemical+r'\b', '<font color="#E65CE6"><b>' + chemical + '</b></font>', pubtator_anns.passages[0].text +" "+ pubtator_anns.passages[1].text, flags=re.I)
    formatted_abstract = re.sub(r'\b'+disease+r'\b', '<font color="#0099FF"><b>' + disease + '</b></font>', formatted_abstract, flags=re.I)


    chemical_from_relation_html = '<font color="#E65CE6"><b>' + chemical + '</b></font>'
    disease_from_relation_html = '<font color="#0099FF"><b>' + disease + '</b></font>'



    # TODO: want something similar to this:
    # paper = relation_task.papers().first()
    # sentences = relation_task.get_sentences()

    ctx = {'chemical' : chemical,
           'disease' : disease,
           'chemical_html': chemical_from_relation_html,
           'disease_html': disease_from_relation_html,
           'current_paper': formatted_abstract,
           'relation': "test"
           }
    return TemplateResponse(request, 'relation/relation.jade', ctx)

form = AnswerForm()

#pass in relation similar to above method
def results(request): #, relation_id
    # relation = Relation.objects.get(pk=relation_id)
    # print request.POST
    #TODO dictionary for max:
    user_work = {'relation': relation,
                 'relation_pair': relation.relation,
                 'relation_type': request.POST['relation_type'],
                 'username': request.user
    }
    # print user_work
    current_answer = Answer.objects.create(relation=relation, relation_pair=relation.relation, relation_type=request.POST['relation_type'], user_confidence='C1', username=request.user)
    relation_type = request.POST['relation_type']

    if request.method == 'POST':
        ctx = { 'form':form,
                'relation_type':relation_type
        }
        return TemplateResponse(request, 'relation/results.jade', ctx)
    #return HttpResponseRedirect(reverse("relation:home"))
