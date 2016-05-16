from __future__ import absolute_import
from django.conf import settings

from configurations import importer
from celery import Celery

import logging
import os
logger = logging.getLogger(__name__)

# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'mark2cure.settings')
os.environ.setdefault('DJANGO_CONFIGURATION', 'Development')

importer.install()

app = Celery('mark2cure')

# Using a string here means the worker will not have to
# pickle the object when using Windows.
app.config_from_object('django.conf:settings')
app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)

app.conf.update(
    CELERY_RESULT_BACKEND='djcelery.backends.database:DatabaseBackend',
)


@app.task(bind=True)
def debug_task(self):
    logger.debug('Debug Task', exc_info=True, extra={'self': self})


