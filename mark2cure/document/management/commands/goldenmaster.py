from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist

from mark2cure.document.models import Document, Section, View, Annotation
import csv, random

class Command(BaseCommand):
    args = '<gm_type>'
    help = 'Import GM documents'

    def handle(self, *args, **options):
        if len(args) < 1: raise Exception('Analysis needs 1 parameters <experiement_id, command>')

        command = args[0]
        self.stdout.write('-- Running GM Routine ({0}) --'.format(command))

        if command == "import":
          self.import_golden_documents()

        elif command == "randomly_make_validation_documents":
          self.randomly_make_validation_documents() #ALREADY RAN ONCE ON PROD

        elif command == "annotate":
          self.annotate_golden_documents()

        else:
          pass

        self.stdout.write('Completed')


    def randomly_make_validation_documents(self):
        documents = Document.objects.filter(source = 'NCBI_corpus_development').all()
        for doc in documents:
          for sec in doc.section_set.all():
            print sec.validate

        # doc_ids = [doc.id for doc in documents]
        # random.shuffle(doc_ids)
        # for doc_id in doc_ids[:10]:
        #     document = Document.objects.get(pk = doc_id)
        #     for section in document.section_set.all():
        #       section.validate = True
        #       section.save()


    def import_golden_documents(self):
        path = "NCBI_corpus_development"
        with open('assets/NCBI_corpus/'+ path +'.txt','r') as f:
            reader = csv.reader(f, delimiter='\t')
            for num, title, text in reader:
                try:
                    doc = Document.objects.get(document_id = num)
                except ObjectDoesNotExist:
                    doc = Document()

                    doc.document_id = num
                    doc.title = title
                    doc.source = path
                    doc.save()

                    sec = Section(kind = "t")
                    sec.text = title
                    sec.document = doc
                    sec.save()

                    sec = Section(kind = "a")
                    sec.text = text
                    sec.document = doc
                    sec.save()


    def annotate_golden_documents(self):
        user, created = User.objects.get_or_create(username="goldenmaster")
        if created:
            user.set_password('')
            user.save()

        path = "NCBI_corpus_development"
        with open('assets/NCBI_corpus/'+ path +'_annos.txt','rU') as f:
            reader = csv.reader(f, delimiter='\t')
            for doc_id, doc_field, ann_type, text, start, stop in reader:
                doc = Document.objects.get(document_id = doc_id)

                for section in doc.section_set.all():
                    if section.kind == doc_field[0]:
                        view, created = View.objects.get_or_create(section = section, user = user)

                        ann, created = Annotation.objects.get_or_create(view = view, text = text, start = start, type = ann_type)
                        ann.kind = "e"
                        ann.user_agent = "goldenmaster"
                        ann.save()


        # Now go back over and confirm they match
        # gm_anns = Annotation.objects.filter(view__user = user).all()
        # for annotation in gm_anns:
        #   text = annotation.view.section.text
        #   print annotation.text, "::", text[annotation.start:]
        #   print "\n - - - - - - \n"




