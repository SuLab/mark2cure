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
import itertools

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
    pubtator_bioc = ""  # need this here for now to return back an empty pubtator_bioc for docs that don't have valid pubtators.

    if current_document.valid_pubtator():
        pubtator_bioc = current_document.as_bioc_with_pubtator_annotations()

        for passage in pubtator_bioc.passages:

            for annotation in passage.annotations:
                concept_type = annotation.infons['type']
                concept_UID = annotation.infons['UID']

                if concept_UID != "None":

                    if concept_type == "0":
                        disease_dict[concept_UID] = annotation.infons

                    if concept_type == "1":
                        gene_dict[concept_UID] = annotation.infons

                    if concept_type == "2":
                        chemical_dict[concept_UID] = annotation.infons

    return disease_dict, gene_dict, chemical_dict, pubtator_bioc


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
            disease_dict, gene_dict, chemical_dict, pubtator_bioc = make_annotation_lists_from_current_document(current_document)
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

    disease_dict, gene_dict, chemical_dict, pubtator_bioc = make_annotation_lists_from_current_document(current_document)

    def make_cgd_annotations(current_document):
        """ This method takes a current document and makes annotations as needed.
        If annotations already exist for this document, then no annotations will be made.
        """
        #TODO this code can be reduced

        if not Annotation.objects.filter(document=current_document).exists():
            if chemical_dict:
                for chemical in chemical_dict:
                    Annotation.objects.create(document=current_document, uid=chemical_dict[chemical]['UID'], stype="c", text=chemical_dict[chemical]['text'], start=0, stop=0)
            if gene_dict:
                for gene in gene_dict:
                    Annotation.objects.create(document=current_document, uid=gene_dict[gene]['UID'], stype="g", text=gene_dict[gene]['text'], start=0, stop=0)
            if disease_dict:
                for disease in disease_dict:
                    Annotation.objects.create(document=current_document, uid=disease_dict[disease]['UID'], stype="d", text=disease_dict[disease]['text'], start=0, stop=0)
            return

        else:
            return

    make_cgd_annotations(current_document)

    def add_relation_pairs_to_database(concept_dict_list, current_document):
        if not Relation.objects.filter(document=current_document).exists():
            # itertools.combinations finds all possible pairs (pairs or 2 here) of relations in any number of lists.
            for tuple_item in list(itertools.combinations(concept_dict_list, 2)):
                concept1_dict = tuple_item[0]
                concept2_dict = tuple_item[1]
                if concept1_dict and concept2_dict:
                    for concept1 in concept1_dict:
                        for concept2 in concept2_dict:
                            concept1_final = Annotation.objects.get(document=current_document, uid=concept1)
                            concept2_final = Annotation.objects.get(document=current_document, uid=concept2)
                            Relation.objects.create(document=current_document, relation=[ concept1, concept2 ], concept1_id=concept1_final, concept2_id=concept2_final, automated_cid=True)
            return
        else:
            return

    concept_dict_list = [gene_dict, chemical_dict, disease_dict]
    add_relation_pairs_to_database(concept_dict_list, current_document)

    # TODO look for annotations instead of relations FASTER TODO TODO

    def find_unanswered_relation(current_document):
        relations = Relation.objects.filter(document=current_document)

        for relation in relations:

            relation_specific_answers = Answer.objects.filter(username=request.user).filter(relation_pair=relation.relation)
            if not relation_specific_answers:
                return relation

    relation = Relation.objects.filter(document=current_document)[0]


    concept1 = Annotation.objects.get(document=current_document, uid=relation.concept1_id)
    concept2 = Annotation.objects.get(document=current_document, uid=relation.concept2_id)

    concept1_text = str(concept1.text)
    concept2_text = str(concept2.text)

    concept1_type = concept1.stype
    concept2_type = concept2.stype

    relation_type = concept1_type + "_" + concept2_type
    # TODO if concept1 and concept2 are stype .... USE different JSON objects for jquery menu....   TODO TODO

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
