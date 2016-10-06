from __future__ import absolute_import


from celery import states
from celery.exceptions import SoftTimeLimitExceeded, Ignore
from ..common import celery_app as app
from celery.schedules import crontab
from celery.task import periodic_task

import logging
logger = logging.getLogger(__name__)


@periodic_task(run_every=crontab(minute='*/5'), ignore_result=True)
def check_system_uptime():
    '''
        Task to run every 5minutes to make sure Celery is running
    '''
    logger.debug('[TASK] check_system_uptime')

