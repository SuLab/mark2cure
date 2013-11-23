from django.db import models
from django.core.exceptions import ObjectDoesNotExist
from django.conf import settings

from mark2cure.document.utils import get_pubmed_document
from pytz import timezone

import datetime

class DocumentManager(models.Manager):
    def pubmed_count(self, pubmed_id):
        return self.filter(document_id__exact = int(pubmed_id)).count()

    def create_from_pubmed_id(self, pubmed_id=None):
        Document = self.model
        pubmed_id = int(pubmed_id)
        print "create_from_pubmed_id"

        ## Check if the account already exists
        try:
            doc = Document.objects.get(document_id = pubmed_id)
            print "found it"
        except ObjectDoesNotExist:
            doc = Document()
            print "had to make it"
            # record = get_pubmed_document(pubmed_id)
            # doc = Document.objects.create(document_id = record.get('PMID'))
            #   doc.title = record.get('TI')
            #   doc.text = record.get('AB')
            #   doc.created = datetime.datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
            #   doc.save()

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
