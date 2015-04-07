from django.dispatch import receiver
from django.db.models.signals import post_save

from .models import Document

import requests


#@receiver(post_save, sender=Document)
#def provider_post_save(sender, instance, created, **kwargs):

def provider_post_save():
    #document = instance
    #if created:
    #    payload = {
    #        'pmid': document.document_id,
    #        'format': 'BioC',
    #        'Disease': 1,
    #        'Gene': 1,
    #        'Chemical': 1,
    #        'Mutation': 0,
    #        'Species': 0}
    #    r = requests.get('http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/PubTator/abstract_ann.cgi', params=payload)

    # Use these in the API:  tmChem,    DNorm,      GNormPlus
    #                       Chemical, Disease,  Gene and Species

    payload = {
        text: 'This is a test',
    }
    requests.get('http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/Chemical/Submit/', params=payload)

    #document.pubtator = r.content
    #document.save()

