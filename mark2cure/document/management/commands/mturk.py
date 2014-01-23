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
        if len(args) < 2: raise Exception('Analysis needs 2 parameters <experiement_id, command>')

        command = args[0]
        documents = args[1]

        self.stdout.write('-- Running MTurk Commands ({0}) on Documents {1} --'.format(command, documents))

        turk = Turk()

        if command == "create_qual":
          turk.make_qualification_test()


        elif command == "disable_all":
          turk.disable_all()


        elif command == "create_hits":
          documents = Document.objects.filter(source = 'NCBI_corpus_development').all()[:10]
          doc_ids = [doc.id for doc in documents]
          random.shuffle(doc_ids)
          # print doc_ids
          for doc_id in doc_ids:
            turk.hit_for_document(doc_id)


        else:
          pass


