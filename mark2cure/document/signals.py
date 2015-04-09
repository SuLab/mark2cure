from django.dispatch import receiver
from django.db.models.signals import post_save

from .models import Document

from mark2cure.common.formatter import bioc_writer
from mark2cure.document.tasks import get_pubtator_response

import requests
from datetime import datetime, timedelta



@receiver(post_save, sender=Document)
def provider_post_save(sender, instance, created, **kwargs):
    doc = instance

    if doc.available_sections().exists() \
            and (doc.pubtator_chem == '' or doc.pubtator_gene == '' or doc.pubtator_disease == ''):

        writer = bioc_writer(None)

        document = doc.as_bioc()
        passage_offset = 0
        for section in doc.available_sections():
            passage = section.as_bioc(passage_offset)
            passage_offset += len(passage.text)
            document.add_passage(passage)
        writer.collection.add_document(document)

        # Use these in the API:  tmChem,    DNorm,      GNormPlus
        #                       Chemical, Disease,  Gene and Species
        payload = {'content-type': 'text/xml'}
        data = str(writer)

        for api_ann in ['tmChem', 'DNorm', 'GNormPlus']:
            response = requests.post(
                    'http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/{api_ann}/Submit/'.format(api_ann=api_ann),
                    data=data, params=payload)

            get_pubtator_response.apply_async(
                args=[
                    doc.pk,
                    api_ann,
                    response.content,
                    data,
                    payload,
                    0],
            )

