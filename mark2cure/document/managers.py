from django.db import models

from mark2cure.common.bioc import BioCReader


class DocumentManager(models.Manager):

    def pubmed_count(self, pubmed_id):
        return self.filter(document_id__exact=int(pubmed_id)).count()

class PubtatorManager(models.Manager):

    def correct_parent_relation(self):
        from mark2cure.document.models import Document
        # Check if each type validates, if so save
        for pubtator in self.all():
            try:
                r = BioCReader(source=pubtator.content)
                r.read()

                pubtator.document = Document.objects.get(document_id=r.collection.documents[0].id)
                pubtator.session_id = ''
            except Exception as e:
                # If one of them doesn't validate leave
                #pubtator.content = None
                print e

            # Do this just so the first time valid_pubtator
            # actually runs we know it's fresh'
            pubtator.validate_cache = False
            pubtator.save()
