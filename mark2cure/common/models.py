from django.contrib.auth.models import User
from django.conf import settings
from django.db import models

from ..document.models import Document, View, Annotation
from ..task.models import DocumentQuestRelationship, Task, UserQuestRelationship

from decimal import Decimal

import random
from collections import Counter

from allauth.account.signals import user_signed_up
from django.dispatch import receiver
from ..task.models import Level
from django.utils import timezone


@receiver(user_signed_up, dispatch_uid='mark2cure.common.allauth.user_signed_up')
def user_signed_up_(request, user, **kwargs):
    task_type_str = request.session.get('initial_training')
    if task_type_str:
        if task_type_str == 'e':
            # After loggin them in, assign the first Level training so we know where to route them
            Level.objects.create(user=user, task_type=request.session.get('initial_training'), level=3, created=timezone.now())
        elif task_type_str == 'r':
            Level.objects.create(user=user, task_type=request.session.get('initial_training'), level=1, created=timezone.now())


class Group(models.Model):
    '''Describe a non-task specific selection of documents (1-n) that curator defined'''

    name = models.CharField(max_length=200)
    stub = models.CharField(max_length=200)
    description = models.TextField(blank=True, null=True)
    order = models.DecimalField(default=0, max_digits=3, decimal_places=3)

    enabled = models.BooleanField(default=False)

    class Meta:
        app_label = 'common'

    def get_documents(self):
        # (TODO?) Return for __in of task_ids
        return Document.objects.filter(task__group=self)

    def total_documents(self):
        # (TODO) rename of return time is reflected
        return DocumentQuestRelationship.objects.filter(task__group=self)

    def doc_count(self):
        dqr_count = len(DocumentQuestRelationship.objects.filter(task__group=self))
        return dqr_count

    def top_five_contributors(self):
        # (TODO) have this return a query set of User models
        """returns a user name list for the group"""
        uqrs = UserQuestRelationship.objects.filter(task__group=self)
        username_list = []
        for uqr in uqrs:
            if uqr.completed:
                # (TODO) not sure if this is the fastest approach
                username_list.append(str.encode(str(uqr.user.username)))
        counter = Counter(username_list)
        top_five_users = [tuple_i[0] for tuple_i in counter.most_common(5)]
        return top_five_users, username_list

    def total_contributors(self):
        """returns a user name list for the group"""
        uqrs = UserQuestRelationship.objects.filter(task__group=self)
        username_list = []
        for uqr in uqrs:
            if uqr.completed:
                # (TODO) not sure if this is the best approach:
                username_list.append(str.encode(str(uqr.user.username)))
        total_contributors = len(set(username_list))
        return total_contributors

    def total_annotation_count(self):
        # (TODO) make work with multiple annotation types
        uqrs = UserQuestRelationship.objects.filter(task__group=self)
        return Annotation.objects.filter(view=View.objects.filter(userquestrelationship=uqrs)).count()

    def current_avg_f(self, weighted=True):
        report_qs = self.report_set.filter(report_type=1).order_by('-created')

        if report_qs.exists():
            report = report_qs.first()
            df = report.dataframe

            if weighted:
                df['wf'] = df['pairings'] * df['f-score']
                try:
                    return df['wf'].sum() / df['pairings'].sum()
                except ZeroDivisionError:
                    return 0.0
            else:
                return df['f-score'].mean()

        else:
            return 0.0

    def percentage_complete(self):
        task_queryset = self.task_set.extra(select={
            "completed": """
                SELECT COUNT(*) AS completed
                FROM task_userquestrelationship
                WHERE (task_userquestrelationship.completed = 1
                    AND task_userquestrelationship.task_id = task_task.id)"""
        })
        completed = sum([x for x in task_queryset.values_list('completed', flat=True) if x is not None])
        required = sum([x for x in task_queryset.values_list('completions', flat=True) if x is not None])
        if required:
            return (Decimal(completed) / Decimal(required)) * 100
        else:
            return 0

    def pubtator_coverage(self):
        """Return back a float representing the amount of documents contained within the ER Group
            have Pubtator annotations

            (TODO) Totally rework this without validate_cache
        """
        # queryset = self.get_documents().prefetch_related('pubtator_set')
        # completed = sum([3 for x in queryset if x.pubtator_set.filter(validate_cache=True).count() == 3])
        # if completed:
        #     return (Decimal(completed) / Decimal(queryset.count() * 3)) * 100
        # else:
        return 0

    def assign(self, documents, smallest_bin=5, largest_bin=5, completions=settings.ENTITY_RECOGNITION_K):
        # Insert the other Documents
        document_set = list(documents.values_list('id', flat=True))
        random.shuffle(document_set)
        last_task = None
        name_counter = 1

        while len(document_set) >= 1:

            quest_size = int(random.uniform(smallest_bin, largest_bin))
            # If there was an existing Task with less than the
            # desired number of documents
            if last_task and last_task.documents.count() < quest_size:

                # Shuffle & Remove the document_pk for use and from being selected again
                random.shuffle(document_set)
                doc_pk = document_set[0]
                document_set.remove(doc_pk)

                document = Document.objects.get(pk=doc_pk)
                DocumentQuestRelationship.objects.create(task=last_task, document=document)

            else:
                last_task, task_created = Task.objects.get_or_create(
                    name=str(name_counter),
                    completions=completions,
                    requires_qualification=7,
                    provides_qualification=7,
                    points=5000,
                    group=self)
                name_counter += 1

    def __unicode__(self):
        return self.name


class SupportMessage(models.Model):
    user = models.ForeignKey(User, blank=True, null=True)
    text = models.TextField()
    referral = models.CharField(max_length=100, blank=True, null=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return u'{text} (via {user})'.format(text=self.text, user=self.user)

    class Meta:
        app_label = 'common'

