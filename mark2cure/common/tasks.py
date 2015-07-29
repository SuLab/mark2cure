from ..common.models import Task, Group, DocumentQuestRelationship
from ..document.models import Document

from celery import task
import random
import logging
logger = logging.getLogger(__name__)


@task()
def check_system_uptime():
    '''
        Task to run every 5minutes to make sure Celery is running
    '''
    logger.info('Celery on rabbitmq-server ran.')


def bin_group(group_pk, document_set_pks):
    group = Group.objects.get(pk=group_pk)

    smallest_bin = 5
    largest_bin = 5
    completions = 15
    random.shuffle(document_set_pks)

    last_task = group.task_set.last()

    while len(document_set_pks) > smallest_bin:

        quest_size = int(random.uniform(smallest_bin, largest_bin))
        # If there was an existing Task with less than the
        # desired number of documents
        if last_task and last_task.documents.count() < quest_size:

            # Shuffle & Remove the document_pk for use and from being selected again
            random.shuffle(document_set_pks)
            doc_pk = document_set_pks[0]
            document_set_pks.remove(doc_pk)

            document = Document.objects.get(pk=doc_pk)
            print 'Add Document', len(document_set_pks), document.valid_pubtator(), last_task
            # if document.valid_pubtator():
            DocumentQuestRelationship.objects.create(task=last_task, document=document)

        else:
            print '> Add New Task'
            if last_task:
                idx = last_task.pk
            else:
                idx = Task.objects.last().pk

            last_task, task_created = Task.objects.get_or_create(
                name=str(idx + 1),
                completions=completions,
                requires_qualification=7,
                provides_qualification=7,
                points=5000,
                group=group)
