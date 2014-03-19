from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
# from django.db.models import Count
from django.contrib.auth.models import User

from mark2cure.document.models import Document, Section, View, Annotation
from mark2cure.document.utils import check_validation_status
from mark2cure.account.models import Ncbo

from mark2cure.common.utils import Turk


import os, os.path, csv, random

class Command(BaseCommand):
    args = '<experiment_run_id>'
    help = 'Command for posting and controlling mturk activity'

    def handle(self, *args, **options):
        if len(args) < 1: raise Exception('Analysis needs 2 parameters <experiement_id, command>')

        command = args[0]
        document_set = ""
        if len(args) > 1:
          document_set = args[1]

        self.stdout.write('-- Running MTurk Commands ({0}) on Documents {1} --'.format(command, document_set))

        turk = Turk()

        if command == "create_qual":
          print turk.make_qualification_test()


        elif command == "disable_all":
          turk.disable_all()


        elif command == "validation":
          # Remove the previous validation rules
          v_sections = Section.objects.filter(validate = True, document__source = 'NCBI_corpus_'+ document_set).all()
          for sec in v_sections:
            sec.validate = False
            sec.save()

          documents = Document.objects.filter(source = 'NCBI_corpus_'+ document_set).all()
          doc_ids = [doc.id for doc in documents]
          random.shuffle(doc_ids)

          for doc_id in doc_ids[:50]:
            secs = Section.objects.filter(document__id = doc_id).all()
            for sec in secs:
              sec.validate = True
              sec.save()


        elif command == "create_hits":
          '''
            Combine all 593 of the training documents with the validation documents from
            the development set (there are 50), randomize their order, then create a HIT
          '''
          documents = Document.objects.filter(source = 'NCBI_corpus_development').all()
          documents = [doc.id for doc in documents]
          random.shuffle(documents)

          print documents
          print len(documents)

          for idx, doc_id in enumerate(documents):
            turk.hit_for_document(doc_id, max_assignments = 10)
            print idx


        else:
          pass


