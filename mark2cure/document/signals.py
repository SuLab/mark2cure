from django.db.models.signals import post_save
from django.dispatch import receiver
from .models import Document, PubtatorRequest

import requests
import re


@receiver(post_save, sender='document.Pubtator')
def pubtator_post_save(sender, instance, created, **kwargs):
    '''
        When a new Pubtator model is created,
        start fetching it's content
    '''
    pubtator = instance

    if created is True and pubtator.content is None and not kwargs.get('raw', False):
        # Make response to post job to pubtator
        payload = {'content-type': 'text/xml'}
        writer = Document.objects.as_writer(documents=[pubtator.document])
        data = str(writer)
        url = 'http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/{api_ann}/Submit/'.format(api_ann=pubtator.kind)

        try:
            response = requests.post(url, data=data, params=payload)
        except Exception as e:
            raise e

        session_id = re.findall(r'\d{4}-\d{4}-\d{4}-\d{4}', response.url)[0]
        PubtatorRequest.objects.get_or_create(
            pubtator=pubtator,
            session_id=session_id)


