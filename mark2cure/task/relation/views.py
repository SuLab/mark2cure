from django.shortcuts import get_object_or_404, redirect
from django.http import HttpResponse
from django.conf import settings
from django.db import connection

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
from .serializers import RelationSerializer, RelationFeedbackSerializer, RelationCereal


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
def relation_task_complete(request, document_pk):
    document = get_object_or_404(Document, pk=document_pk)
    first_section = document.section_set.first()
    view = get_object_or_404(View, task_type='ri', completed=True, section=first_section, user=request.user)

    ctx = {
        'document': document,
        'points': request.user.profile.score(view=view)
    }
    return TemplateResponse(request, 'relation/task-complete.jade', ctx)


@login_required
def submit_document_set(request, document_pk):
    document = get_object_or_404(Document, pk=document_pk)

    if request.method == 'POST':

        first_section = document.section_set.first()
        view = View.objects.get(task_type='ri', completed=False, section=first_section, user=request.user)
        view.completed = True
        view.save()

        view_content_type = ContentType.objects.get_for_model(view)

        Point.objects.create(user=request.user,
                             amount=settings.RELATION_DOC_POINTS,
                             content_type=view_content_type,
                             object_id=view.id,
                             created=timezone.now())

        return redirect('task-relation:task-complete', document_pk=document.pk)


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
                             amount=settings.RELATION_REL_POINTS,
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
def document_analysis(request, document_pk, relation_pk=None):
    """
        API for returning analysis details for Document
        Relation task completions
    """
    document = get_object_or_404(Document, pk=document_pk)

    relation = None
    # If a relation was specified, only show results for that
    if relation_pk:
        relation = get_object_or_404(Relation, pk=relation_pk)

    # Start the DB Connection
    c = connection.cursor()

    cmd_str = """
        SELECT  `relation_relation`.`id` as `relationship_id`,
                `relation_relation`.`document_id`,
                `relation_relation`.`relation_type`,
                `relation_relation`.`concept_1_id`,
                `relation_relation`.`concept_2_id`,

                # ANY_VALUE(`concept_relationship_1`.`stype`) as `concept_1_stype`,
                ANY_VALUE(`concept_text_1`.`text`) as `concept_1_text`,

                # ANY_VALUE(`concept_relationship_2`.`stype`) as `concept_2_stype`,
                ANY_VALUE(`concept_text_2`.`text`) as `concept_2_text`

        FROM `relation_relation`

        INNER JOIN `relation_conceptdocumentrelationship` as `concept_relationship_1`
            ON `concept_relationship_1`.`document_id` = `relation_relation`.`document_id`
            INNER JOIN `relation_concepttext` as `concept_text_1`
                    ON `concept_text_1`.`concept_id` = `relation_relation`.`concept_1_id`

        INNER JOIN `relation_conceptdocumentrelationship` as `concept_relationship_2`
            ON `concept_relationship_2`.`document_id` = `relation_relation`.`document_id`
            INNER JOIN `relation_concepttext` as `concept_text_2`
                    ON `concept_text_2`.`concept_id` = `relation_relation`.`concept_2_id`

        WHERE `relation_relation`.`document_id` = {document_id} {relation_logic}
        GROUP BY `relation_relation`.`id`
    """.format(
        document_id=document.pk,
        relation_logic=' AND `relation_relation`.`id` = {0}'.format(relation.pk) if relation else '')

    c.execute(cmd_str)
    rel_tasks = []
    for x in c.fetchall():
        rel_tasks.append(x)

    cmd_str = """
        SELECT  `document_document`.`id` as `doc_pk`,
                `document_document`.`document_id` as `pmid`,
                `relation_relationannotation`.`relation_id` as `relationship_id`,
                `relation_relationannotation`.`answer` as `relationship_answer`,
                `document_annotation`.`created`,
                `document_view`.`user_id`

        FROM `relation_relationannotation`

        INNER JOIN `document_annotation`
            ON `document_annotation`.`object_id` = `relation_relationannotation`.`id` AND `document_annotation`.`kind` = 'r'

            INNER JOIN `document_view`
                ON `document_annotation`.`view_id` = `document_view`.`id`

                INNER JOIN `document_section`
                    ON `document_section`.`id` = `document_view`.`section_id`

                INNER JOIN `document_document`
                    ON `document_document`.`id` = `document_section`.`document_id`

        WHERE `document_document`.`id` =  {document_id}
    """.format(document_id=document.pk)
    c.execute(cmd_str)

    rel_submissions = []
    for x in c.fetchall():
        rel_submissions.append(x)

    from collections import defaultdict
    groups = defaultdict(list)
    for obj in rel_submissions:
        groups[obj[2]].append(obj)

    serializer = RelationCereal(rel_tasks, many=True, context={'sub_dict': groups, 'user': request.user})
    return Response(serializer.data)


@login_required
@api_view(['GET'])
def fetch_relation_feedback(request, relation_pk):
    relation = get_object_or_404(Relation, pk=relation_pk)
    serializer = RelationFeedbackSerializer([relation], many=True, context={'user': request.user})
    return Response(serializer.data[0])

