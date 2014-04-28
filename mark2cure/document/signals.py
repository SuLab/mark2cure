from django.conf import settings
from django.db.models import signals

from mark2cure.document.models import Document

import requests, importlib, logging
logger = logging.getLogger(__name__)



def document_save_handler(sender, instance, **kwargs):
    '''
      Automatically seed new documents onto MTurk
    '''
    document = instance
    pass


signals.post_save.connect(document_save_handler, sender=Document)


