from django.dispatch import receiver
from django.db.models.signals import post_save

from .models import Document, Pubtator

from mark2cure.document.tasks import get_pubtator_response

from datetime import datetime, timedelta
import requests


def fetchPubtator(sender, instance, created, **kwargs):
    '''
        When a new Pubtator model is created,
        start fetching it's content
    '''
    pubtator = instance

    if created:
        # Make response to post job to pubtator
        payload = {'content-type': 'text/xml'}
        writer = pubtator.as_writer()
        data = str(writer)
        url = 'http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/{api_ann}/Submit/'.format(api_ann=pubtator.kind)
        response = requests.post(url, data=data, params=payload)
        pubtator.session_id = response.content
        pubtator.save()

        # (TODO) This is bad but pk passing wasn't working
        get_pubtator_response.apply_async(
            args=[pubtator, data, payload, 0],
        )

post_save.connect(fetchPubtator, sender=Pubtator)


@receiver(post_save, sender=Document)
def provider_post_save(sender, instance, created, **kwargs):
    '''
      When a new document is made, add Bioc documents for
      the types we want to collect
    '''
    doc = instance

    if created:
        for api_ann in ['tmChem', 'DNorm', 'GNormPlus']:
            #if doc.available_sections().exists() and not Pubtator.objects.filter(document=doc, kind=api_ann).exists():
            pub, pub_c = Pubtator.objects.get_or_create(document=doc, kind=api_ann)
            pub.save(force_update=True)
