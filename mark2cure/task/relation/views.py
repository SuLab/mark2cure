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
from .serializers import DocumentRelationSerializer, RelationAnalysisSerializer


@login_required
def relation_task(request, document_pk):
        """Document base page for completing available Relations.

            1) Fetches the Document Relations API to get available work for the document
            2) Displays the Tree selection interface for comparing the 2 concepts
                - Submits currents and shows next over AJAX
            3) When relationships are completed, flag the associated View as completed
                - This flag is not used to determine anything
        """
        document = get_object_or_404(Document, pk=document_pk)

        if request.POST:
            """API to trigger completion of the Document for the user to the extent we allow

                We award them Points for completing the available Relationships
                that are currently available to them. They may recieve this bonus
                multiple times

                (TODO) Block out a user from multiple completions of Document
            """
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

        # (TODO) Return if the user has already completed 20 Relations within this Document
        ctx = {
            'document': document,
        }
        return TemplateResponse(request, 'relation/task.jade', ctx)


@login_required
def relation_task_complete(request, document_pk):
    """ Document conclusion page for thanking b/c the user exhausted the total number
        of relation submissions we allow (20)

        - Allows entrance to Talk Page or Return to Dashboard
    """
    document = get_object_or_404(Document, pk=document_pk)
    first_section = document.section_set.first()
    view = get_object_or_404(View, task_type='ri', completed=True, section=first_section, user=request.user)

    # (TODO) This should be changed from 'r' to concept type
    anns_user_doc_list = Annotation.objects.filter(kind='r', view_id__user_id=request.user.id,
                                                   view_id__section__document_id=document_pk)

    reln_id_list = []
    relation_list = []
    reln_ann_list = []
    for ann in anns_user_doc_list:

        reln_ann = RelationAnnotation.objects.get(id=ann.object_id)
        reln = Relation.objects.get(id=reln_ann.relation_id)
        if reln.id not in reln_id_list:
            relation_list.append(reln)
            reln_ann_list.append(reln_ann)
        reln_id_list.append(reln.id)

    ctx = {
        'document_pk': document.pk,
        'reln_ann_list': reln_ann_list,
        'document': document,
        'points': request.user.profile.score(view=view)
    }
    return TemplateResponse(request, 'relation/task-complete.jade', ctx)


@login_required
@require_http_methods(['POST'])
def submit_annotation(request, document_pk, relation_pk):
    """ Submit the selected relation type as an Annotation that is associated
        to a View of the first Document Section (.first() if Task Type isn't Section Specific)
    """
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
    """ API that includes all the Relations within for a Document for the user (TODO)
    """
    document = get_object_or_404(Document, pk=document_pk)

    # Create the views
    for section in document.section_set.all():
        View.objects.get_or_create(task_type='ri', section=section, user=request.user)

    cmd_str = ""
    with open('mark2cure/task/relation/commands/get-user-relations-for-document.sql', 'r') as f:
        cmd_str = f.read()
    cmd_str = cmd_str.format(
        document_id=document.pk,
        user_id=request.user.pk,
        completions=settings.ENTITY_RECOGNITION_K)

    # Start the DB Connection
    c = connection.cursor()
    c.execute(cmd_str)

    queryset = [{'id': x[0],
                 'document_id': x[1],
                 'relation_type': x[2],
                 'progress': x[3],
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

    serializer = DocumentRelationSerializer(queryset, many=True, context={'user': request.user})
    return Response(serializer.data)


@login_required
@api_view(['GET'])
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
    for x in groups.iterkeys():
        item = groups[x][0]
        item['answers'] = groups[x]
        queryset.append(item)

    serializer = RelationAnalysisSerializer(queryset, many=True)
    return Response(serializer.data)

