from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType
from django.conf import settings

from rest_framework.generics import ListCreateAPIView
from rest_framework.response import Response
from rest_framework import status, permissions

from django.shortcuts import get_object_or_404, redirect
from django.views.decorators.http import require_http_methods
from django.http import HttpResponse, HttpResponseServerError
from django.template.response import TemplateResponse

from ...common.formatter import bioc_as_json, apply_annotations, clean_df
from ...document.models import Document, Annotation
from ..models import Level, Task, UserQuestRelationship
from .models import EntityRecognitionAnnotation
from .utils import generate_results, select_best_opponent
from ...score.models import Point
from .serializers import AnnotationSerializer

from django.utils import timezone

import random


@login_required
def user_pmid_results_bioc(request, doc_pk, user_pk, format_type):
    '''
        Return a BioC file for the PMID with only a specific user's M2C entity_recognition submissions as annotations
    '''
    document = get_object_or_404(Document, pk=doc_pk)
    user = get_object_or_404(User, pk=user_pk)

    df = Document.objects.entity_recognition_df(documents=[document], users=[user], include_pubtator=False)
    df = clean_df(df)

    # BioC Writer Response that will serve all partner comparison information
    writer = Document.objects.as_writer(documents=[document])
    writer = apply_annotations(writer, er_df=df)

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')


@login_required
@require_http_methods(['POST'])
def identify_annotations_results_bioc(request, task_pk, doc_pk, format_type):
    '''
        Returns back the modified BioC file with additional collection metadata that allow us to compare you to a user and show score data
            (TODO) Expand this section
    '''

    task = get_object_or_404(Task, pk=task_pk)
    document = task.documents.filter(pk=doc_pk).first()
    if not document:
        return HttpResponseServerError()

    '''
        Try to find an optimal opponete to pair the player
        against. If one isn't available or none meet the minimum
        requirements then just tell the player they've
        annotated a new document
    '''
    opponent = select_best_opponent(task, document, request.user)
    if not opponent:
        # No other work has ever been done on this apparently
        # so we reward the user and let them know they were
        # first via a different template / bonus points
        uqr = task.userquestrelationship_set.filter(user=request.user).first()

        # Did the new user provide at least 1 annotation?
        # (TODO) Did the new annotations differ from pubtator?
        # it would make more sense to do that as a valid check of
        # "contribution" effort

        content_type = ContentType.objects.get_for_model(task)
        Point.objects.create(user=request.user,
                             amount=settings.ENTITY_RECOGNITION_DOC_POINTS,
                             content_type=content_type,
                             object_id=task.id,
                             created=timezone.now())
        return HttpResponseServerError('points_awarded')

    df = Document.objects.entity_recognition_df(documents=[document], users=[opponent], include_pubtator=False)
    df = clean_df(df)

    # BioC Writer Response that will serve all partner comparison information
    writer = Document.objects.as_writer(documents=[document])
    writer = apply_annotations(writer, er_df=df)

    # Other results exist if other people have at least viewed
    # the quest and we know other users have at least submitted
    # results for this particular document
    player_views = []
    opponent_views = []
    for section in document.available_sections():
        # If paired against a player who has completed the task multiple times
        # compare the to the first instance of the person completing that Document <==> Quest
        # while taking the latest version of the player's

        uqr = task.userquestrelationship_set.filter(user=request.user).first()
        player_view = uqr.views.filter(user=request.user, section=section, completed=True).first()

        quest_rel = task.userquestrelationship_set.filter(user=opponent).first()
        opponent_view = quest_rel.views.filter(section=section, completed=True).first()

        player_views.append(player_view)
        opponent_views.append(opponent_view)

        # Save who the player was paired against
        player_view.opponent = opponent_view
        player_view.save()

    results = generate_results(player_views, opponent_views)
    score = results[0][2] * settings.ENTITY_RECOGNITION_DOC_POINTS
    if score > 0:
        Point.objects.create(
            user=request.user,
            amount=score,
            content_type=ContentType.objects.get_for_model(task),
            object_id=task.id,
            created=timezone.now())

    writer.collection.put_infon('flatter', random.choice(settings.POSTIVE_FLATTER) if score > 500 else random.choice(settings.SUPPORT_FLATTER))
    writer.collection.put_infon('points', str(int(round(score))))
    writer.collection.put_infon('partner', str(opponent.username))
    writer.collection.put_infon('partner_level', Level.objects.filter(user=opponent, task_type='e').first().get_name())

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')


@login_required
def quest_read_doc(request, quest_pk, doc_idx):
    task = get_object_or_404(Task, pk=quest_pk)

    # Redirect if trying to access more documents
    # than this task contains
    if int(doc_idx) > task.remaining_documents_count():
        return redirect('task-entity-recognition:quest-home', quest_pk=task.pk)

    # Confirm the user has started, but not completed this Task
    user_quest_relationship = task.user_relationship(request.user, False)
    if not user_quest_relationship:
        return redirect('task-entity-recognition:quest-home', quest_pk=task.pk)

    task_doc_pks_completed = user_quest_relationship.completed_document_ids()
    if int(doc_idx) <= len(task_doc_pks_completed):
        return redirect('task-entity-recognition:quest-home', quest_pk=task.pk)

    '''
    # only take the one that has views
    user_quest_rel = task.userquestrelationship_set.filter(user=user, completed=False).latest()
    user_quest_rel_views = user_quest_rel.views

    # If a completed view doesn't exist, redirect them to complete the document submission
    if not user_quest_rel_views.filter(section__document=doc, completed=True).exists():
        return redirect('document:read', task.pk, doc.pk)
    '''

    # Fetch available documents
    task_doc_uncompleted = task.remaining_documents(task_doc_pks_completed)
    random.shuffle(task_doc_uncompleted)

    document = task_doc_uncompleted[0]
    # task.create_views(document, request.user)

    ctx = {'task': task,
           'completed_doc_pks': task_doc_pks_completed,
           'uncompleted_docs': task_doc_uncompleted,
           'document': document}
    return TemplateResponse(request, 'entity_recognition/quest.jade', ctx)


@login_required
def quest_read_doc_results_bioc(request, quest_pk, doc_pk, user_pk, format_type):
    '''
        Returns the BioC file for how a particular user did
        on a document WITHIN a specific quest
    '''
    # Check pass just for url santization reasons
    get_object_or_404(Task, pk=quest_pk)
    document = get_object_or_404(Document, pk=doc_pk)
    user = get_object_or_404(User, pk=user_pk)

    df = Document.objects.entity_recognition_df(documents=[document], users=[user], include_pubtator=False)
    df = clean_df(df)

    # BioC Writer Response that will serve all partner comparison information
    writer = Document.objects.as_writer(documents=[document])
    writer = apply_annotations(writer, er_df=df)

    writer.collection.put_infon('partner', str(user.username))
    writer.collection.put_infon('partner_level', Level.objects.filter(user=user, task_type='e').first().get_name())

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')


@login_required
def quest_read_doc_results(request, quest_pk, doc_idx):
    '''
        Allows player to revist result page to look at comparision
    '''
    task = get_object_or_404(Task, pk=quest_pk)
    # Get the UQR. Completed=False b/c they may want to view results
    # of previous doc_idx's before finishing the rest of the document
    user_quest_relationship = task.user_relationship(request.user, True)

    # The completed document at this index
    # Using Abstracts but this [1,1,2,2,3,3,4,4,5,5] assumption
    # is extremely fragile
    relevant_views = user_quest_relationship.views.filter(section__kind='a')

    if relevant_views.filter(opponent__isnull=False).exists():
        relevant_view = relevant_views.filter(opponent__isnull=False)[int(doc_idx)]
        opponent = relevant_view.opponent.user
    else:
        # Novel annotations for the player (unpaired)
        relevant_view = relevant_views.all()[int(doc_idx)]
        opponent = None

    ctx = {'task': task,
           'opponent': opponent,
           'document': relevant_view.section.document}
    return TemplateResponse(request, 'entity_recognition/quest-results.jade', ctx)


class SubmitDocumentApiView(ListCreateAPIView):
    """API to submit annotations and mark the document as completed
    """
    # request, quest_pk, document_pk
    serializer_class = AnnotationSerializer(many=True)
    permission_classes = (permissions.IsAuthenticated,)

    def create(self, request, *args, **kwargs):
        data = request.data
        task = get_object_or_404(Task, pk=self.kwargs['quest_pk'])
        document = get_object_or_404(Document, pk=self.kwargs['document_pk'])

        serializer = AnnotationSerializer(data=data, many=True)
        if serializer.is_valid():
            user_quest_rel = task.user_relationship(request.user, False)

            if not user_quest_rel:
                return HttpResponseServerError('User Quest Relationship not found')

            for d in data:
                view = user_quest_rel.views.filter(section_id=d.get('section_pk'), completed=False).first()

                if not view:
                    return HttpResponseServerError('View for Annotation not found')

                er_ann = EntityRecognitionAnnotation.objects.create(
                    type_idx=d.get('type_id'),
                    text=d.get('text'),
                    start=d.get('start')
                )
                er_ann_content_type = ContentType.objects.get_for_model(er_ann)
                Annotation.objects.create(
                    kind='e',
                    view=view,
                    content_type=er_ann_content_type,
                    object_id=er_ann.pk)

            for view in user_quest_rel.views.filter(section__document=document):
                view.completed = True
                view.save()

            headers = self.get_success_headers(serializer.data)
            return Response(serializer.data, status=status.HTTP_201_CREATED, headers=headers)


@login_required
def quest_submit(request, task, bypass_post=False):
    # (TODO) Add validation check here at some point

    if request.POST or bypass_post:
        user_quest_relationship = task.user_relationship(request.user, False)

        if not user_quest_relationship.completed:
            Point.objects.create(
                user=request.user,
                amount=task.points,
                content_type=ContentType.objects.get_for_model(task),
                object_id=task.id,
                created=timezone.now())

        user_quest_relationship.completed = True
        user_quest_relationship.save()


@login_required
def quest_feedback(request, quest_pk):
    task = get_object_or_404(Task, pk=quest_pk)
    ctx = {'task': task}
    return TemplateResponse(request, 'entity_recognition/quest-feedback.jade', ctx)


@login_required
def quest_read(request, quest_pk):
    task = get_object_or_404(Task, pk=quest_pk)

    # Check if user has pre-existing relationship with Quest
    user_quest_rel_queryset = UserQuestRelationship.objects.filter(task=task, user=request.user)

    if user_quest_rel_queryset.exists():

        # User has viewed this Quest before and Completed it
        # so show them the feedback page
        if user_quest_rel_queryset.filter(completed=True).exists():
            return redirect('task-entity-recognition:quest-feedback', quest_pk=task.pk)

        # Not using get_or_create b/c get occasionally returned multiple (unknown bug source)
        user_quest_relationship = user_quest_rel_queryset.first()
        task_doc_pks_completed = user_quest_relationship.completed_document_ids()

        # If there are no more documents to do, mark the Quest
        # as completed and go to dashboard
        task_doc_uncompleted = task.remaining_documents(task_doc_pks_completed)
        if len(task_doc_uncompleted) == 0:
            quest_submit(request, task, True)
            return redirect('task-entity-recognition:quest-feedback', quest_pk=task.pk)

        next_doc_idx = len(task_doc_pks_completed) + 1
        return redirect('task-entity-recognition:quest-document', quest_pk=task.pk, doc_idx=next_doc_idx)

    else:
        # Create the User >> Quest relationship
        UserQuestRelationship.objects.create(task=task, user=request.user, completed=False)

        documents = list(task.documents.all())
        for document in documents:
            task.create_views(document, request.user)

        # Route the user to the right idx doc
        return redirect('task-entity-recognition:quest-document', quest_pk=task.pk, doc_idx=1)

