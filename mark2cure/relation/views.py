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
from ..document.models import Document, Section, Pubtator

from brabeion import badges
import datetime
import re
import json
import random

from .models import Answer, Relation
from .forms import AnswerForm
import os


#TODO, needs improvement
def make_annotation_lists_from_current_document(current_document):
    """
    TODO: should be removed from views and possibly added to models or tasks?
    Input is one document. Output is three concept lists in order: disease, gene,
    chemical. Also output is the pubtator_anns.

    This is SO slow*** Need to find a way to only show things if there are chemicals + diseases.

    If there are no chem/dis, then no reason to show it.
    """
    disease_dict = {}
    gene_dict = {}
    chemical_dict = {}
    pubtator_anns = ''
    if current_document.valid_pubtator():
        pubtator_anns = current_document.as_bioc_with_pubtator_annotations_jf()
        num_passages = len(pubtator_anns.passages)
        for i in range(0, num_passages):
            len_annotations = len(pubtator_anns.passages[i].annotations)
            for k in range(0, len_annotations):
                concept_type = pubtator_anns.passages[i].annotations[k].infons['type']
                concept_UID = pubtator_anns.passages[i].annotations[k].infons['UID']
                if concept_UID != "None":
                    if concept_type == "0":
                        disease_dict[concept_UID] = pubtator_anns.passages[i].annotations[k].infons
                    if concept_type == "1":
                        gene_dict[concept_UID] = pubtator_anns.passages[i].annotations[k].infons
                    if concept_type == "2":
                        chemical_dict[concept_UID] = pubtator_anns.passages[i].annotations[k].infons
    return disease_dict, gene_dict, chemical_dict, pubtator_anns


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
                user_paper_answers = Answer.objects.all().filter(relation_pair=j.relation).filter(username=request.user)
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
            disease_dict, gene_dict, chemical_dict, pubtator_anns = make_annotation_lists_from_current_document(current_document)
            num_non_empty_dicts = 0
            if not len(disease_dict) == 0:
                num_non_empty_dicts += 1
            if not len(gene_dict) == 0:
                num_non_empty_dicts += 1
            if not len(chemical_dict) == 0:
                num_non_empty_dicts += 1
            if num_non_empty_dicts < 2:
                exclude_list.append(i)

        return exclude_list

    exclude_list = exclude_empty_concept_lists(queryset_papers) + exclude_user_answered_relations(queryset_papers)

    exclude_these_list = [exclude_this.id for exclude_this in exclude_list]
    filtered_queryset_papers = Document.objects.all().exclude(id__in=exclude_these_list)[:30]

    ctx = {
        'documents': filtered_queryset_papers,
    }
    return TemplateResponse(request, 'relation/home.jade', ctx)


@login_required
def relation(request, document_pk):
    form = AnswerForm
    current_document = get_object_or_404(Document, pk=document_pk)

    #relations = Relation.objects.filter(document=current_document.id)

    disease_dict, gene_dict, chemical_dict, pubtator_anns = make_annotation_lists_from_current_document(current_document)
    # for relation in relations:
    #     # relations that user has already completed
    #     relation_specific_answers = Answer.objects.filter(username=request.user).filter(relation_pair=relation.relation)
    #     if not relation_specific_answers:
    #         break

    try:
        chemical_UID = chemical_dict.keys()[0]
    except:
        pass
    try:
        disease_UID = disease_dict.keys()[0]
    except:
        pass
    try:
        gene_UID = gene_dict.keys()[0]
    except:
        pass

    #print disease_UID, gene_UID, chemical_UID, "hello"

    def make_cgd_annotations(current_document):
        for chemical in chemical_dict:
            Annotation.objects.create(document=current_document, uid=chemical_dict[chemical]['UID'], stype="c", text=chemical_dict[chemical]['text'], start=0, stop=0)
        for gene in gene_dict:
            Annotation.objects.create(document=current_document, uid=gene_dict[gene]['UID'], stype="g", text=gene_dict[gene]['text'], start=0, stop=0)
        for disease in disease_dict:
            Annotation.objects.create(document=current_document, uid=disease_dict[disease]['UID'], stype="d", text=disease_dict[disease]['text'], start=0, stop=0)

    make_cgd_annotations(current_document)

    relation = Relation.objects.create(document=current_document, relation=[ chemical_dict[chemical_UID]['UID'], disease_dict[disease_UID]['UID'] ], chemical_id=chemical_dict[chemical_UID]['UID'], disease_id=disease_dict[disease_UID]['UID'], automated_cid=True)
    relation = Relation.objects.filter(document=current_document)[0]


    chemical = chemical_dict[chemical_UID]['text']
    disease = disease_dict[disease_UID]['text']


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
           'current_document': current_document,
           'relation': relation
           }
    return TemplateResponse(request, 'relation/relation.jade', ctx)


#pass in relation similar to above method
def results(request, relation_id): #, relation_id
    relation = get_object_or_404(Relation, pk=relation_id)

    relation_type = request.POST['relation_type']

    form = AnswerForm(request.POST or None)

    # print form
    if form.is_valid():
        save_it = form.save(commit=False)
        save_it = save_it.save()

    if request.method == 'POST':
        ctx = {
        'relation_type':relation_type
        }
        return TemplateResponse(request, 'relation/results.jade', ctx)
    #return HttpResponseRedirect(reverse("relation:home"))
