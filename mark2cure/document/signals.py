from django.db.models.signals import post_save
from django.dispatch import receiver

from .models import Pubtator

import requests


@receiver(post_save, sender=Pubtator)
def pubtator_post_save(sender, instance, created, **kwargs):
    '''
        When a new Pubtator model is created,
        start fetching it's content
    '''
    pubtator = instance

    if created and pubtator.content is None:
        # Make response to post job to pubtator
        payload = {'content-type': 'text/xml'}
        writer = pubtator.document.as_writer()
        data = str(writer)
        url = 'http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/{api_ann}/Submit/'.format(api_ann=pubtator.kind)

        response = requests.post(url, data=data, params=payload)
        pubtator.session_id = response.content
        pubtator.save()

