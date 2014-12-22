from django.db import models


class DocumentManager(models.Manager):

    def pubmed_count(self, pubmed_id):
        return self.filter(document_id__exact=int(pubmed_id)).count()

