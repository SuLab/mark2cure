from django.contrib.contenttypes.models import ContentType
from django.contrib.auth.decorators import login_required
from django.contrib.contenttypes.models import ContentType
from django.template.response import TemplateResponse
from django.shortcuts import get_object_or_404

from django_comments.models import Comment

from .decorators import doc_completion_required
from ..document.models import Document, Annotation
from ..task.entity_recognition.models import EntityRecognitionAnnotation

from django.contrib import messages
from django.utils import timezone
from collections import Counter
from datetime import timedelta


@login_required
@doc_completion_required
def home(request, pubmed_id):
    document = get_object_or_404(Document, document_id=pubmed_id)
    content_type_id = str(ContentType.objects.get_for_model(EntityRecognitionAnnotation.objects.all().first()).id)
    disease = Counter(EntityRecognitionAnnotation.objects.annotations_texts_for_document_and_type(document.pk, 'disease', content_type_id))
    gene_protein = Counter(EntityRecognitionAnnotation.objects.annotations_texts_for_document_and_type(document.pk, 'gene_protein', content_type_id))
    drug = Counter(EntityRecognitionAnnotation.objects.annotations_texts_for_document_and_type(document.pk, 'drug', content_type_id))

    ctx = {
        'doc': document,
        'diseases': disease.most_common(20),
        'gene_proteins': gene_protein.most_common(20),
        'drugs': drug.most_common(20)
    }

    return TemplateResponse(request, 'talk/home.jade', ctx)


@login_required
def annotation_search(request):
    annotation = request.GET.get('q')
    completed_document_pks = request.user.profile.completed_document_pks()

    # (TODO) Sanitize search param for RAW SQL
    content_type_id = str(ContentType.objects.get_for_model(EntityRecognitionAnnotation.objects.all().first()).id)
    document_pks = EntityRecognitionAnnotation.objects.document_pks_by_text_and_document_pks(annotation, completed_document_pks, content_type_id)

    documents = [(group[1], Document.objects.get(pk=group[0])) for group in Counter(document_pks).most_common()]
    ctx = {
        'annotation': annotation,
        'documents': documents
    }
    return TemplateResponse(request, 'talk/annotation_search.jade', ctx)


@login_required
def recent_discussion(request):
    doc_content_pk = ContentType.objects.get(name='document').pk
    completed_document_pks = request.user.profile.completed_document_pks()

    is_moderator = request.user.groups.filter(name='Comment Moderators').exists()
    last_week = timezone.now().date() - timedelta(days=7)

    content_type_id = str(ContentType.objects.get_for_model(EntityRecognitionAnnotation.objects.all().first()).id)
    if is_moderator:
        msg = '<p class="lead text-center">You\'re a moderator, showing Global View.</p>'
        messages.info(request, msg, extra_tags='safe alert-info')

        comment_queryset = Comment.objects.filter(
            content_type_id=doc_content_pk)

        annotations = EntityRecognitionAnnotation.objects.annotations_texts_by_created(last_week, content_type_id)

        comment_documents_pks = Comment.objects.filter(
            content_type_id=doc_content_pk).values_list('object_pk', flat=True)

    else:
        comment_queryset = Comment.objects.filter(
            content_type_id=doc_content_pk,
            object_pk__in=completed_document_pks)

        annotations = EntityRecognitionAnnotation.objects.annotations_texts_by_created_and_document_pks(last_week, completed_document_pks, content_type_id)

        comment_documents_pks = Comment.objects.filter(
            content_type_id=doc_content_pk,
            object_pk__in=completed_document_pks).values_list('object_pk', flat=True)

    recent_comments = comment_queryset.extra(select={
        "pmid": """
            SELECT document_document.document_id AS pmid
            FROM document_document
            WHERE document_document.id = django_comments.object_pk
            LIMIT 1"""
    }).order_by('-submit_date')

    annotation_counter = Counter(annotations)

    documents = Document.objects.filter(
        pk__in=comment_documents_pks
    ).extra(select={
        "comment_count": """
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
        'documents': documents,
    }
    return TemplateResponse(request, 'talk/recent_discussion.jade', ctx)

