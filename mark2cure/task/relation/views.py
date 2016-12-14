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

from .models import Relation, RelationAnnotation, ConceptDocumentRelationship
from .serializers import RelationSerializer, RelationCereal


class RelationTask(View):

    @login_required
    # @method_decorator(user_passes_test(lambda u: u.is_superuser))
    # @method_decorator(user_passes_test(lambda u: u.is_staff))
    def dispatch(self, *args, **kwargs):
        return super(RelationTask, self).dispatch(*args, **kwargs)

    def get(self, request, *args, **kwargs):
        """Document base page for completing available Relations.
        """
        document = get_object_or_404(Document, pk=kwargs['document_pk'])

        # (TODO) Return if the user has already completed 20 Relations within this Document

        ctx = {
            'document': document,
        }
        return TemplateResponse(request, 'relation/task.jade', ctx)

    def post(self, request, *args, **kwargs):
        """API to trigger completion of the Document for the user to the extent we allow
        """
        document = get_object_or_404(Document, pk=kwargs['document_pk'])

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
    """ API for returning analysis details for Document
        Relation task completions
    """
    document = get_object_or_404(Document, pk=document_pk)

    relation = None
    # If a relation was specified, only show results for that
    if relation_pk:
        relation = get_object_or_404(Relation, pk=relation_pk)

    cmd_str = ""
    with open('mark2cure/task/relation/commands/get-relations-for-document.sql', 'r') as f:
        cmd_str = f.read()
    cmd_str = cmd_str.format(
        document_id=document.pk,
        relation_logic=' AND `relation_relation`.`id` = {0}'.format(relation.pk) if relation else '')

    # Start the DB Connection
    c = connection.cursor()
    c.execute(cmd_str)

    rel_tasks = []

    # (TODO) Ugly patch for #189 (https://github.com/SuLab/mark2cure/issues/189)
    # Intended to mirror behavior at mark2cure/task/relation/serializers.pyL55
    cdr_query = ConceptDocumentRelationship.objects.filter(document=document)
    for x in c.fetchall():
        cdr1 = cdr_query.filter(concept_text__concept_id=x[3]).first()
        cdr2 = cdr_query.filter(concept_text__concept_id=x[4]).first()
        rel_tasks.append((x[0], x[1], x[2], x[3], x[4], cdr1.concept_text.text, cdr2.concept_text.text))

    cmd_str = ""
    with open('mark2cure/task/relation/commands/get-relations-annotations-for-document.sql', 'r') as f:
        cmd_str = f.read()
    cmd_str = cmd_str.format(document_id=document.pk)

    c.execute(cmd_str)

    rel_submissions = []
    for x in c.fetchall():
        rel_submissions.append(x)

    # Close the connection
    c.close()

    from collections import defaultdict
    groups = defaultdict(list)
    for obj in rel_submissions:
        groups[obj[2]].append(obj)

    serializer = RelationCereal(rel_tasks, many=True, context={'sub_dict': groups, 'user': request.user})
    return Response(serializer.data)

