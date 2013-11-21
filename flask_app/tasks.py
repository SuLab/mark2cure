# -*- coding: utf-8 -*-
"""
    mark2cure.tasks
    ~~~~~~~~~~~~~~

    mark2cure tasks module
"""

from .core import mail
from .factory import create_celery_app

celery = create_celery_app()


@celery.task
def send_manager_added_email(*recipients):
    print 'sending manager added email...'


@celery.task
def send_manager_removed_email(*recipients):
    print 'sending manager removed email...'

