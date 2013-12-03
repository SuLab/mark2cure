from django.db import models
from django.core.exceptions import ObjectDoesNotExist
from django.conf import settings

# from mark2cure.document.models import Section

from pytz import timezone
from Bio import Entrez, Medline

import datetime

class DocumentManager(models.Manager):
    def pubmed_count(self, pubmed_id):
        return self.filter(document_id__exact = int(pubmed_id)).count()

    def create_from_pubmed_id(self, pubmed_id=None):
        Document = self.model
        pubmed_id = str(pubmed_id)

        ## Check if the account already exists
        try:
            doc = Document.objects.get(document_id = pubmed_id)
        except ObjectDoesNotExist:
            doc = Document()

            Entrez.email = settings.ENTREZ_EMAIL
            h = Entrez.efetch(db='pubmed', id=[pubmed_id], rettype='medline', retmode='text')
            records = Medline.parse(h)

            for record in records:
              if record.get('TI') and record.get('AB') and record.get('PMID') and record.get('CRDT'):
                doc.document_id = record.get('PMID')
                doc.title = record.get('TI')
                doc.created = datetime.datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
                doc.source = "pubmed"
                doc.save()

                # sec = Section()

                # doc.title = record.get('TI')
                # doc.text = record.get('AB')
              break

        return doc


    def get_random_document(self):
        '''
        This is documented as being potentially expensive, we may want to do something like
        http://stackoverflow.com/a/6405601 instead
        '''
        return self.order_by('?')[0]


class AnnotationManager(models.Manager):
    pass
#
#     def clean_annotations(self):
#         Annotation = self.model
