from django.contrib.auth.decorators import login_required
from django.contrib.contenttypes.models import ContentType
from django.template.response import TemplateResponse
from django.shortcuts import get_object_or_404

from django_comments.models import Comment

from .decorators import doc_completion_required
from ..document.models import Document
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
    disease = Counter(EntityRecognitionAnnotation.objects.annotations_texts_for_document_and_type(document.pk, 0, content_type_id))
    gene_protein = Counter(EntityRecognitionAnnotation.objects.annotations_texts_for_document_and_type(document.pk, 1, content_type_id))
    drug = Counter(EntityRecognitionAnnotation.objects.annotations_texts_for_document_and_type(document.pk, 2, content_type_id))

    ctx = {
        'doc': document,
        'diseases': disease.most_common(20),
        'gene_proteins': gene_protein.most_common(20),
        'drugs': drug.most_common(20)
    }

    return TemplateResponse(request, 'talk/home.html', ctx)


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
    return TemplateResponse(request, 'talk/annotation_search.html', ctx)


@login_required
def recent_discussion(request):
    doc_content_pk = ContentType.objects.get(app_label='document', model='document').pk
    completed_document_pks = request.user.profile.completed_document_pks()

    is_moderator = request.user.groups.filter(name='Comment Moderators').exists()
    last_week = timezone.now().date() - timedelta(days=7)

    content_type_id = str(ContentType.objects.get_for_model(EntityRecognitionAnnotation.objects.all().first()).id)


    if is_moderator:
        msg = '<p class="lead text-xs-center">You\'re a moderator, showing Global View.</p>'
        messages.info(request, msg, extra_tags='safe alert-info')


    ctx = {
        'comments': recent_comments,
        'documents': documents,
    }
    return TemplateResponse(request, 'talk/recent_discussion.html', ctx)

