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
from ..userprofile.models import UserProfile
from ..document.models import Document, Section, Pubtator

from brabeion import badges
import datetime
import re
import json
import random

from .models import Answer, Relation, Concept
from .forms import AnswerForm
import os


#TODO, needs improvement



@login_required
def home(request):
    # how many documents to display on the home page, sort through 30 docs and see if there are docs that have chem/diseases
    queryset_papers = Document.objects.all()[:30]

    # exclude if user already answered
    def exclude_user_answered_relations(queryset_papers):
        exclude_list = []
        for i in queryset_papers:
            user_paper_relations = Relation.objects.filter(document=i.id)
            flag = 0
            for j in user_paper_relations:
                user_paper_answers = Answer.objects.all().filter(relation=j.pk).filter(username=request.user)
                if len(user_paper_answers) == 1:
                    flag += 1
                if flag == len(user_paper_relations):
                    exclude_list.append(i)

        return exclude_list # returns a query_set

    ### exclude if 2 or 3 lists are empty
    def exclude_empty_concept_lists(queryset_papers):
        exclude_list = []
        for i in queryset_papers:
            current_document = Document.objects.get(pk=i.pk)
            disease_dict, gene_dict, chemical_dict, pubtator_bioc = current_document.make_concept_lists()
            total_full_dicts = 0
            total_full_dicts = len(disease_dict) + len(gene_dict) + len(chemical_dict)
            if total_full_dicts < 2:
                exclude_list.append(i)

        return exclude_list

    exclude_list = exclude_empty_concept_lists(queryset_papers) + exclude_user_answered_relations(queryset_papers)

    exclude_these_list = [exclude_this.id for exclude_this in exclude_list]
    filtered_queryset_papers = Document.objects.all().exclude(id__in=exclude_these_list)[:10]

    ctx = {
        'documents': filtered_queryset_papers,
    }
    return TemplateResponse(request, 'relation/home.jade', ctx)


@login_required
def relation(request, document_pk):
    form = AnswerForm
    current_document = get_object_or_404(Document, pk=document_pk)

    disease_dict, gene_dict, chemical_dict, pubtator_bioc = current_document.make_concept_lists()
    concept_dict_list = [gene_dict, chemical_dict, disease_dict]
    current_document.make_cgd_concepts(disease_dict, gene_dict, chemical_dict)

    current_document.add_relation_pairs_to_database(concept_dict_list)

    relation = current_document.unanswered_relation(request)

    unanswered_relations_for_user = current_document.unanswered_relation_list(request)

    concept1 = Concept.objects.get(document=current_document, uid=relation.concept1_id)
    concept2 = Concept.objects.get(document=current_document, uid=relation.concept2_id)

    concept1_text = str(concept1.text)
    concept2_text = str(concept2.text)

    concept1_type = concept1.stype
    concept2_type = concept2.stype

    # TODO add this to the model (relation.stype)
    relation_type = concept1_type + "_" + concept2_type
    # TODO if concept1 and concept2 are stype .... USE different JSON objects for jquery menu....   TODO TODO


    def make_relation_dict(unanswered_relations_for_user):
        relation_list = []

        for relation in unanswered_relations_for_user:

            concept1 = Concept.objects.get(document=current_document, uid=relation.concept1_id)
            concept2 = Concept.objects.get(document=current_document, uid=relation.concept2_id)

            concept1_text = str(concept1.text)
            concept2_text = str(concept2.text)

            concept1_type = concept1.stype
            concept2_type = concept2.stype

            relation_type = concept1_type + "_" + concept2_type


            relation_list.append({relation: {'concept1': concept1,
                                       'concept2': concept2,
                                       'concept1_text': concept1_text,
                                       'concept2_text': concept2_text,
                                       'concept1_type': concept1_type,
                                       'concept2_type': concept2_type,
                                       'relation_type': relation_type} })
        return relation_list


    relation_list = make_relation_dict(unanswered_relations_for_user)

    formatted_abstract = re.sub(r'\b' + concept1_text + r'\b', '<font color="#E65CE6"><b>' + concept1_text + '</b></font>', pubtator_bioc.passages[0].text +" "+ pubtator_bioc.passages[1].text, flags=re.I)
    formatted_abstract = re.sub(r'\b' + concept2_text + r'\b', '<font color="#0099FF"><b>' + concept2_text + '</b></font>', formatted_abstract, flags=re.I)

    chemical_from_relation_html = '<font color="#E65CE6"><b>' + concept1_text + '</b></font>'
    disease_from_relation_html = '<font color="#0099FF"><b>' + concept2_text + '</b></font>'

    # TODO: want something similar to this:
    # paper = relation_task.papers().first()
    # sentences = relation_task.get_sentences()

    ctx = {'concept1' : concept1_text,
           'concept2' : concept2_text,
           'chemical_html': chemical_from_relation_html,
           'disease_html': disease_from_relation_html,
           'current_paper': formatted_abstract,
           'current_document': current_document,
           'relation_type': relation_type,
           'relation': relation,
           'relation_list': relation_list
           }
    return TemplateResponse(request, 'relation/relation.jade', ctx)

#pass in relation similar to above method
def results(request, relation_id): #, relation_id
    relation = get_object_or_404(Relation, pk=relation_id)
    relation_type = request.POST['relation_type']
    form = AnswerForm(request.POST or None)

    if form.is_valid():
        form.save()

    if request.method == 'POST':
        ctx = {
        'relation_type':relation_type
        }
        return TemplateResponse(request, 'relation/results.jade', ctx)
    #return HttpResponseRedirect(reverse("relation:home"))

@login_required
@require_http_methods(['POST'])
def create_post(request):
    form = AnswerForm(request.POST or None)
    if form.is_valid():
        form.save()

    return HttpResponse(200)
