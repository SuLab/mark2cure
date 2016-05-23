from django.shortcuts import get_object_or_404
from django.http import HttpResponse
from django.conf import settings

from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.template.response import TemplateResponse
from django.contrib.contenttypes.models import ContentType
from django.utils import timezone

from rest_framework.decorators import api_view
from rest_framework.response import Response

from ...score.models import Point
from ...document.models import Document, View, Annotation

from .models import Relation, RelationAnnotation
from .serializers import RelationSerializer, RelationFeedbackSerializer


@login_required
def relation_task_home(request, document_pk):
    """
        Main page for users to find the relationship between two concepts.
    """
    document = get_object_or_404(Document, pk=document_pk)

    ctx = {
        'document': document,
    }
    return TemplateResponse(request, 'relation/task.jade', ctx)


@login_required
def show_document_results(request, document_pk):
    document = get_object_or_404(Document, pk=document_pk)

    if request.method == 'POST':

        first_section = document.section_set.first()
        view = View.objects.get(task_type='ri', completed=False, section=first_section, user=request.user)
        view.completed = True
        view.save()

        ctx = {
            'document': document,
        }
        return TemplateResponse(request, 'relation/results.jade', ctx)


@login_required
@require_http_methods(['POST'])
def submit_annotation(request, document_pk, relation_pk):
    '''
        Submit the selected relation type as an Annotation that is associated
        to a View of the first Document Section (.first() if Task Type isn't Section Specific)
    '''
    document = get_object_or_404(Document, pk=document_pk)
    relation = get_object_or_404(Relation, pk=relation_pk)

    first_section = document.section_set.first()
    view = View.objects.get(task_type='ri', completed=False, section=first_section, user=request.user)

    current_selection = request.POST.get('relation', None)
    if current_selection:
        relation_ann = RelationAnnotation.objects.create(relation=relation, answer=current_selection)
        relation_ann_content_type = ContentType.objects.get_for_model(relation_ann)

        Annotation.objects.create(
            kind='r',
            view=view,
            content_type=relation_ann_content_type,
            object_id=relation_ann.id)

        # Assign a point to the specific Relation Annotation
        Point.objects.create(user=request.user,
                             amount=settings.RELATION_DOC_POINTS,
                             content_type=relation_ann_content_type,
                             object_id=relation_ann.id,
                             created=timezone.now())

        return HttpResponse(200)

    return HttpResponse(500)


@login_required
@api_view(['GET'])
def fetch_document_relations(request, document_pk):
    """API that includes all the relationships within one document based on the
    document's primary key.
    """
    document = get_object_or_404(Document, pk=document_pk)

    # Create the views
    for section in document.section_set.all():
        View.objects.get_or_create(task_type='ri', section=section, user=request.user)

    queryset = Relation.objects.filter(document=document).extra(select={
        "user_completed_count": """
            SELECT COUNT(*) AS user_completed_count
            FROM document_annotation

            INNER JOIN `document_view` ON ( document_annotation.view_id = document_view.id )
            INNER JOIN `relation_relationannotation` ON ( document_annotation.object_id = relation_relationannotation.id )

            WHERE (document_annotation.kind = 'r'
                AND document_annotation.content_type_id = 56
                AND document_view.user_id = %d
                AND relation_relationannotation.relation_id = relation_relation.id)""" % (request.user.pk,)
    })

    serializer = RelationSerializer(queryset, many=True, context={'user': request.user})
    return Response(serializer.data)


@login_required
@api_view(['GET'])
def fetch_relation_feedback(request, relation_pk):
    """
        API for returning total number of answers for a relation.
    """
    relation = get_object_or_404(Relation, pk=relation_pk)

    serializer = RelationFeedbackSerializer([relation], many=True, context={'user': request.user})

    return Response(serializer.data[0])
