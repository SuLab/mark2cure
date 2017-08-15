from django.shortcuts import get_object_or_404
from django.http import HttpResponse
from django.conf import settings
from django.db import connection

from django.contrib.auth.decorators import login_required
from django.template.response import TemplateResponse
from django.contrib.contenttypes.models import ContentType
from django.utils import timezone

from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status

from ...score.models import Point
from ...document.models import Document, View, Annotation

from .models import Relation, RelationAnnotation
from .serializers import DocumentRelationSerializer, RelationAnalysisSerializer


@login_required
@api_view(['GET', 'POST'])
def document_analysis(request, document_pk, relation_pk=None):
    """ API for returning analysis details for Document
        Relation task completions

        Scoped only towards relationships the user is aware of
    """
    document = get_object_or_404(Document, pk=document_pk)

    relation_logic = ''
    # If a relation was specified, only show results for that
    if relation_pk:
        relation = get_object_or_404(Relation, pk=relation_pk)
        relation_logic = 'AND `relation_relation`.`id` = {0}'.format(relation.pk)

    cmd_str = ""
    with open('mark2cure/task/relation/commands/get-relations-analysis-for-user-and-document.sql', 'r') as f:
        cmd_str = f.read()
    cmd_str = cmd_str.format(
        document_id=document.pk,
        user_id=request.user.pk,
        relation_logic=relation_logic
    )

    # Start the DB Connection
    c = connection.cursor()
    c.execute(cmd_str)

    queryset = [{
        'id': x[0],
        'document_id': x[1],
        'kind': x[2],
        'concept_1_id': x[3],
        'concept_1_text': x[4],

        'concept_2_id': x[5],
        'concept_2_text': x[6],

        'answer_hash': x[7],
        'user_id': x[8],
        'self': x[9],
    } for x in c.fetchall()]

    # Close the connection
    c.close()

    from collections import defaultdict
    groups = defaultdict(list)
    for obj in queryset:
        groups[obj.get('id')].append(obj)

    queryset = []
    for x in groups.keys():
        item = groups[x][0]
        item['answers'] = groups[x]
        queryset.append(item)

    serializer = RelationAnalysisSerializer(queryset, many=True)
    return Response(serializer.data)


@login_required
@api_view(['GET'])
def re_task_relationships_list(request, document_pk):
    """ API that includes all the Relations within for a Document for the user (TODO)
    """
    document = get_object_or_404(Document, pk=document_pk)

    # Create the views
    for section in document.section_set.all():
        View.objects.get_or_create(task_type='ri', section=section, user=request.user)

    cmd_str = ""
    with open('mark2cure/task/relation/commands/get-user-relations-for-document.sql', 'r') as f:
        cmd_str = f.read()

    # Start the DB Connection
    c = connection.cursor()

    c.execute('SET @k_max = {completions};'.format(completions=settings.ENTITY_RECOGNITION_K))
    c.execute('SET @user_id = {user_id};'.format(user_id=request.user.pk))
    c.execute('SET @document_id = {document_id};'.format(document_id=document.pk))
    c.execute(cmd_str)

    queryset = [{'id': x[0],
                 'document_id': x[1],
                 'relation_type': x[2],

                 'community_progress': x[3],
                 'community_completed': x[4],
                 'user_completed': x[5],

                 'concept_1_id': x[6],
                 'concept_1_type': x[7],
                 'concept_1_text': x[8],

                 'concept_2_id': x[9],
                 'concept_2_type': x[10],
                 'concept_2_text': x[11]} for x in c.fetchall()]

    # Close the connection
    c.close()

    serializer = DocumentRelationSerializer(queryset, many=True)
    return Response(serializer.data)


@login_required
@api_view(['POST'])
def re_task_relationship_submit(request, document_pk, relation_pk):
    """ Submit the selected relation type as an Annotation that is associated
        to a View of the first Document Section (.first() if Task Type isn't Section Specific)
    """
    document = get_object_or_404(Document, pk=document_pk)
    relation = get_object_or_404(Relation, pk=relation_pk)

    first_section = document.section_set.first()
    view = View.objects.get(task_type='ri', completed=False, section=first_section, user=request.user)

    current_selection = request.POST.get('relation', None)
    if current_selection:
        cmd_str = ""
        with open('mark2cure/task/relation/commands/get-relation-ann-exists.sql', 'r') as f:
            cmd_str = f.read()
        cmd_str = cmd_str.format(
            document_id=document.pk,
            user_id=request.user.pk,
            relation_id=relation.pk
        )

        # Start the DB Connection
        c = connection.cursor()
        c.execute(cmd_str)
        already_submitted = c.fetchone()[0]

        if already_submitted:
            content = {'message': 'A Relationship Extraction Annotation can only be submitted once.'}
            return Response(content, status=status.HTTP_409_CONFLICT)

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

        return document_analysis(request, document_pk, relation_pk)

    return HttpResponse(500)


@login_required
@api_view(['POST'])
def re_task_submit(request, document_pk):
    '''Method to POST to in order to confirm all Document RE work
        was completed for a Task and the RE Task should now be considered complete.

        returns REDocumentResult (tree.js)
        -- -- -- --
        API to trigger completion of the Document for the user to the extent we allow
        We award them Points for completing the available Relationships
        that are currently available to them. They may recieve this bonus
        multiple times
        (TODO) Block out a user from multiple completions of Document
    '''
    document = get_object_or_404(Document, pk=document_pk)
    first_section = document.section_set.first()
    view = View.objects.filter(task_type='ri', section=first_section, user=request.user).last()

    if view.completed:
        re_task_created = False
        award = Point.objects.filter(user=request.user,
                                     content_type=ContentType.objects.get_for_model(View),
                                     object_id=view.id).first()

    else:
        re_task_created = True
        award = Point.objects.create(user=request.user,
                                     amount=settings.RELATION_DOC_POINTS,
                                     content_type=ContentType.objects.get_for_model(view),
                                     object_id=view.id)
        view.completed = True
        view.save()

    return Response({
        'document': {
            'pk': document.pk,
            'pmid': document.document_id,
            'title': document.title,
            'relationship_count': Relation.objects.filter(document=document).count()
        },

        're_task': {  # The View
            'pk': view.pk,
            'created': re_task_created
        },

        'award': {
            'pk': award.pk,
            'amount': int(round(award.amount))
        },

    })


@login_required
def re_task(request, document_pk):
    '''View that serves required HTML and starts the Tree library

        Document shuffling, and other RE document ordering logic
        is done by the client.
    '''
    document = get_object_or_404(Document, pk=document_pk)

    # (TODO) Return if the user has already completed 20 Relations within this Document

    return TemplateResponse(request, 'relation/task.jade', {'document': document})


