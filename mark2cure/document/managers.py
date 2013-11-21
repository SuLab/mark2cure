from django.db import models
from django.core.exceptions import ObjectDoesNotExist
from django.conf import settings
from pytz import timezone
import datetime

class DocumentManager(models.Manager):

    def create_from_pubmed_id(self, pubmed_id=None):
        Document = self.model

        ##########
        #
        # (TODO) Do MEDLINE magic here
        #
        ##########

        ## Check if the account already exists
        try:
            doc = Document.objects.get(pubmed_id = pubmed_id)
        except ObjectDoesNotExist:
            doc = Document()

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
