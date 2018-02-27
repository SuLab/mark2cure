from django.contrib.auth.decorators import login_required
from django.contrib.contenttypes.models import ContentType
from django.template.response import TemplateResponse
from django.shortcuts import get_object_or_404

from .decorators import doc_completion_required
from ..document.models import Document
from ..task.ner.models import EntityRecognitionAnnotation

from collections import Counter


@login_required
@doc_completion_required
def home(request, document_pk):
    document = get_object_or_404(Document, pk=document_pk)
    return TemplateResponse(request, 'talk/home.html', {'document': document})


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
    return TemplateResponse(request, 'talk/recent_discussion.html')

