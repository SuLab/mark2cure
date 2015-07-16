from django.template.response import TemplateResponse
from django.shortcuts import get_object_or_404, redirect
from django.conf import settings
from django.views.decorators.http import require_http_methods

from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.contrib.messages import get_messages
from django.contrib import messages
from django.http import HttpResponse

from django.contrib.contenttypes.models import ContentType
from django_comments.models import Comment

from mark2cure.userprofile.models import UserProfile
from mark2cure.document.models import Document, Annotation

from .decorators import doc_completion_required

from collections import Counter

from datetime import timedelta
from django.utils import timezone

import logging
import random
import os
logger = logging.getLogger(__name__)


@login_required
def annotation_search(request):
    annotation = request.GET.get('q')
    completed_document_pks = request.user.profile.completed_document_pks()
    queryset = Annotation.objects.filter(
            text=annotation,
            view__section__document__pk__in=completed_document_pks
        ).values_list('view__section__document__pk', flat=True)

    documents = [(group[1], Document.objects.get(pk=group[0])) for group in Counter(queryset).most_common()]
    ctx = {
        'annotation': annotation,
        'documents': documents
    }
    return TemplateResponse(request, 'talk/annotation_search.jade', ctx)

@login_required
def recent_discussion(request):
    doc_content_pk = ContentType.objects.get(name='document').pk
    completed_document_pks = request.user.profile.completed_document_pks()

    recent_comments = Comment.objects.filter(content_type_id=doc_content_pk, object_pk__in=completed_document_pks).extra(select = {
        "pmid" : """
            SELECT document_document.document_id AS pmid
            FROM document_document
            WHERE document_document.id = django_comments.object_pk
            LIMIT 1"""
    }).order_by('-submit_date')

    last_week = timezone.now().date() - timedelta(days=7)
    annotations = Annotation.objects.filter(
            created__gte=last_week,
            view__section__document__pk__in=completed_document_pks
        ).exclude(text='').values_list('text', flat=True)
    annotation_counter = Counter(annotations)

    comment_documents_pks = Comment.objects.filter(content_type_id=doc_content_pk, object_pk__in=completed_document_pks).values_list('object_pk', flat=True)
    documents = Document.objects.filter(
            pk__in=comment_documents_pks
        ).extra(select = {
        "comment_count" : """
            SELECT count(*) AS comment_count
            FROM django_comments
            WHERE (django_comments.content_type_id = 23
                AND django_comments.object_pk = document_document.id)"""
    })
    documents = list(documents)
    documents.sort(key=lambda x: x.comment_count, reverse=True)

    ctx = {
        'comments': recent_comments,
        'annotations': annotation_counter.most_common(100),
        'documents': documents
    }
    return TemplateResponse(request, 'talk/recent_discussion.jade', ctx)


@login_required
@doc_completion_required
def home(request, pubmed_id):
    document = get_object_or_404(Document, document_id=pubmed_id)

    disease = Counter( document.annotations().filter(type='disease').exclude(text='').values_list('text', flat=True) )
    gene_protein = Counter( document.annotations().filter(type='gene_protein').exclude(text='').values_list('text', flat=True) )
    drug = Counter( document.annotations().filter(type='drug').exclude(text='').values_list('text', flat=True) )

    ctx = {
        'doc': document,
        'diseases': disease.most_common(20),
        'gene_proteins': gene_protein.most_common(20),
        'drugs': drug.most_common(20)
    }

    return TemplateResponse(request, 'talk/home.jade', ctx)


