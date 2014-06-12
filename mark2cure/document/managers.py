from django.db import models


class DocumentManager(models.Manager):


    def pubmed_count(self, pubmed_id):
        return self.filter(document_id__exact = int(pubmed_id)).count()


    def get_random_document(self):
        '''
        This is documented as being potentially expensive, we may want to do something like
        http://stackoverflow.com/a/6405601 instead
        '''
        return self.order_by('?')[0]


