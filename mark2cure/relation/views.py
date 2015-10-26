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
import itertools

from .models import Answer, Relation, Concept
from .forms import AnswerForm
import os


@login_required
def home(request):
    # how many documents to display on the home page, sort through 30 docs and see if there are docs that have chem/diseases
    # TODO test
    concepts_available = Concept.objects.all()

    queryset_documents = []
    for concept in concepts_available:
        if concept.document not in queryset_documents:
            queryset_documents.append(concept.document)

    # TODO exclude documents that do not contain unanswered user relation pairs
    def exclude_user_answered_relations(queryset_documents):
        exclude_list = []
        for document in queryset_documents:
            user_paper_relations = Relation.objects.filter(document=document.id)
            flag = 0
            for relation in user_paper_relations:
                user_paper_answers = Answer.objects.all().filter(relation=relation.pk).filter(username=request.user)
                if len(user_paper_answers) == 1:
                    flag += 1
                if flag == len(user_paper_relations):
                    exclude_list.append(document)

        return exclude_list # returns a query_set

    exclude_list = exclude_user_answered_relations(queryset_documents)

    for doc in exclude_list:
        if doc in queryset_documents:
            queryset_documents.remove(doc)

    # Just show a few
    queryset_documents = queryset_documents[:10]

    ctx = {
        'documents': queryset_documents,
    }
    return TemplateResponse(request, 'relation/home.jade', ctx)


@login_required
def relation(request, document_pk):
    form = AnswerForm
    # the document instance
    current_document = get_object_or_404(Document, pk=document_pk)
    # the section instances
    section_title = Section.objects.get(document=current_document, kind="t")
    section_abstract = Section.objects.get(document=current_document, kind="a")

    relation = current_document.unanswered_relation(request)
    unanswered_relations_for_user = current_document.unanswered_relation_list(request)

    concept1 = Concept.objects.get(document=current_document, uid=relation.concept1_id)
    concept2 = Concept.objects.get(document=current_document, uid=relation.concept2_id)

    relation_type = concept1.stype + "_" + concept2.stype

    def make_relation_dict(unanswered_relations_for_user):
        relation_list = []

        for relation in unanswered_relations_for_user:

            concept1 = Concept.objects.get(document=current_document, uid=relation.concept1_id)
            concept2 = Concept.objects.get(document=current_document, uid=relation.concept2_id)

            relation_type = concept1.stype + "_" + concept2.stype

            relation_list.append({relation.pk: { 'concept1': concept1,
                                                 'concept2': concept2,
                                                 'relation_type': relation_type} })
        return relation_list

    relation_list = make_relation_dict(unanswered_relations_for_user)

    # TODO make my own api for this view

    pk_list = []
    for dict_item in relation_list:
        pk_list.append(dict_item.keys())

    pk_list = list(itertools.chain(*pk_list))
    pk_list = json.dumps(pk_list)

    formatted_abstract = re.sub(r'\b' + str(concept1.text) + r'\b', '<font color="#E65CE6"><b>' + str(concept1.text) + '</b></font>', str(section_title) +" "+ str(section_abstract), flags=re.I)
    formatted_abstract = re.sub(r'\b' + str(concept2.text) + r'\b', '<font color="#0099FF"><b>' + str(concept2.text) + '</b></font>', formatted_abstract, flags=re.I)

    chemical_from_relation_html = '<font color="#E65CE6"><b>' + str(concept1.text) + '</b></font>'
    disease_from_relation_html = '<font color="#0099FF"><b>' + str(concept2.text) + '</b></font>'

    # TODO: want something similar to this:
    # paper = relation_task.papers().first()
    # sentences = relation_task.get_sentences()

    print concept1
    print concept2
    print chemical_from_relation_html
    print disease_from_relation_html
    print relation_type
    print relation_list

    ctx = {'concept1' : concept1,
           'concept2' : concept2,
           'chemical_html': chemical_from_relation_html,
           'disease_html': disease_from_relation_html,
           'current_paper': formatted_abstract,
           'current_document': current_document,
           'relation_type': relation_type,
           'relation': relation,
           'relation_list': relation_list,
           'pk_list': pk_list
           }
    return TemplateResponse(request, 'relation/relation3.jade', ctx)

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

@login_required
@require_http_methods(['POST'])
def jen_bioc(request):
    # print request
    # print "HELLO"
    relation_list = request.POST['relation_list']
    # print relation_list
    writer_json = bioc_as_json(relation_list)
    return HttpResponse(200)
