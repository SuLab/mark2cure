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
    logger.debug('[TASK] check_system_uptime')

