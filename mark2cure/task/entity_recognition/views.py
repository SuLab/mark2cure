from django.contrib.auth.models import User
from django.contrib.auth.decorators import login_required
from django.contrib.contenttypes.models import ContentType
from django.conf import settings

from rest_framework.generics import ListCreateAPIView
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status, permissions

from django.shortcuts import get_object_or_404
from django.http import HttpResponseServerError
from django.template.response import TemplateResponse

from ...document.models import Document, Annotation
from ...common.models import Group
from ..models import Level, Task, UserQuestRelationship
from .models import EntityRecognitionAnnotation
from .utils import generate_results, select_best_opponent
from ...score.models import Point
from .serializers import AnnotationSerializer

from django.utils import timezone

import random


@login_required
@api_view(['GET'])
def ner_quest_document_results(request, task_pk, doc_pk):
    """Returns back submission results that allow a comparison to an opponent
    """
    task = get_object_or_404(Task, pk=task_pk)
    document = task.documents.filter(pk=doc_pk).first()
    if not document:
        return HttpResponseServerError()

    opponent_pk = select_best_opponent(task_pk, doc_pk, request.user.pk)
    opponent_anns = {}

    if opponent_pk:
        opponent_user = User.objects.get(pk=opponent_pk)
        opponent_dict = {
            'pk': opponent_pk,
            'name': opponent_user.username,
            'level': Level.objects.filter(user_id=opponent_pk, task_type='e').first().get_name(),
        }

        player_view_pks = []
        opponent_view_pks = []
        for section in document.available_sections():
            uqr = task.userquestrelationship_set.filter(user=request.user).first()
            player_view = uqr.views.filter(user=request.user, section=section, completed=True).first()
            player_view_pks.append(player_view.pk)

            quest_rel = task.userquestrelationship_set.filter(user=opponent_user).first()
            opponent_view = quest_rel.views.filter(section=section, completed=True).first()
            opponent_view_pks.append(opponent_view.pk)

            # Save who the player was paired against
            player_view.opponent = opponent_view
            player_view.save()

        results = generate_results(player_view_pks, opponent_view_pks)
        points = results[0][2] * settings.ENTITY_RECOGNITION_DOC_POINTS  # F Score * Point Multiplier (1000)

        opponent_ner_ann_df = Document.objects.ner_df(document_pks=[doc_pk], user_pks=[opponent_pk], include_pubtator=False)
        for section_id, group in opponent_ner_ann_df.groupby('section_id'):
            opponent_anns[int(section_id)] = []
            for ann_group_idx, ann in group.iterrows():
                opponent_anns[int(section_id)].append({
                    'type_id': int(ann['ann_type_idx']),
                    'start': int(ann['start_position']),
                    'text': str(ann['text'])
                })

    else:
        opponent_dict = None
        # No other work has ever been done on this apparently
        # so we reward the user and let them know they were
        # first via a different template / bonus points
        uqr = task.userquestrelationship_set.filter(user=request.user).first()
        points = settings.ENTITY_RECOGNITION_DOC_POINTS

    # print('Completed?', uqr.views.filter(section__document_id=doc_pk).values_list('completed', flat=True) )
    # if user_quest_relationship.completed:
    #     uqr_created = False
    #     award = Point.objects.filter(
    #         content_type=ContentType.objects.get_for_model(task),
    #         object_id=task.id,
    #         user=request.user).first()
    # else:
    #     uqr_created = True
    #     award = Point.objects.create(
    #         user=request.user,
    #         amount=task.points,
    #         content_type=ContentType.objects.get_for_model(task),
    #         object_id=task.id,
    #         created=timezone.now())

    award = Point.objects.create(
        user=request.user,
        amount=points,
        content_type=ContentType.objects.get_for_model(task),
        object_id=task.id,
        created=timezone.now())

    res = {
        'task_pk': task.pk,
        'document_pk': document.pk,
        'flatter': random.choice(settings.POSTIVE_FLATTER) if award.amount > 500 else random.choice(settings.SUPPORT_FLATTER),
        'award': {
            'pk': award.pk,
            'amount': int(round(award.amount))
        },
        'opponent': opponent_dict,
        'opponent_anns': opponent_anns
    }

    return Response(res)


class NERDocumentSubmissionView(ListCreateAPIView):
    '''API to submit annotations and mark the document as completed
    '''
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
@api_view(['POST'])
def ner_quest_submit(request, quest_pk):
    '''Method to POST to in order to confirm all Document NER work
        was completed for a Quest and the Quest should now be considered complete.
    '''
    # (TODO) Add validation check here at some point
    # (TODO) Return a dict for a NERQuestComplete model

    task = get_object_or_404(Task, pk=quest_pk)
    group = get_object_or_404(Group, pk=task.group_id)
    user_quest_relationship = UserQuestRelationship.objects.filter(task=task, user=request.user).last()

    if user_quest_relationship.completed:
        uqr_created = False
        award = Point.objects.filter(
            content_type=ContentType.objects.get_for_model(task),
            object_id=task.id,
            user=request.user).first()

    else:
        uqr_created = True
        award = Point.objects.create(
            user=request.user,
            amount=task.points,
            content_type=ContentType.objects.get_for_model(task),
            object_id=task.id,
            created=timezone.now())

        user_quest_relationship.completed = True
        user_quest_relationship.save()

    return Response({
        'task': {
            'pk': task.pk,
            'name': task.name
        },
        'group': {
            'pk': group.pk,
            'name': group.name,
            'stub': group.stub,
            'enabled': group.enabled,
            'description': group.description
        },
        'uqr_pk': user_quest_relationship.pk,
        'uqr_created': uqr_created,
        'award': {
            'pk': award.pk,
            'amount': int(round(award.amount))
        }
    })


@login_required
@api_view(['GET'])
def ner_quest(request, quest_pk):
    '''View that serves required HTML and starts the YPet library

        Document shuffling, and other NER document ordering logic
        is done by the client.
    '''
    task = get_object_or_404(Task, pk=quest_pk)

    # Check if user has pre-existing relationship with Quest
    if not UserQuestRelationship.objects.filter(task=task, user=request.user).exists():
        # Create the User >> Quest relationship
        UserQuestRelationship.objects.create(task=task, user=request.user, completed=False)

        documents = list(task.documents.all())
        for document in documents:
            task.create_views(document, request.user)

    return TemplateResponse(request, 'entity_recognition/quest.jade', {'task_pk': task.pk})
