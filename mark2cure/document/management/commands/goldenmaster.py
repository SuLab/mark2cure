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
        if len(args) < 2: raise Exception('Analysis needs 1 parameters <experiement_id, command>')

        command = args[0]
        document_set = args[1]
        self.stdout.write('-- Running GM Routine ({0}) --'.format(command))

        if command == "import":
          self.import_golden_documents(document_set)

        elif command == "randomly_make_validation_documents":
          self.randomly_make_validation_documents(document_set) #ALREADY RAN ONCE ON PROD

        elif command == "annotate":
          self.annotate_golden_documents(document_set)

        else:
          pass

        self.stdout.write('Completed')




    def import_golden_documents(self, document_set):

        with open('assets/datasets/'+ document_set +'_cleaned.txt','r') as f:
            reader = csv.reader(f, delimiter='\t')
            for num, title, text in reader:
                print title

                doc, doc_c = Document.objects.get_or_create(document_id = num)
                doc.title = title
                doc.source = document_set
                doc.save()

                sec, sec_c = Section.objects.get_or_create(kind = "t", document = doc)
                sec.text = title
                sec.save()

                sec, sec_c = Section.objects.get_or_create(kind = "a", document = doc)
                sec.text = text
                sec.save()



    def randomly_make_validation_documents(self, document_set):
        documents = Document.objects.filter(source = document_set).all()
        # for doc in documents:
        #   for sec in doc.section_set.all():
        #     print sec.validate

        # doc_ids = [doc.id for doc in documents]
        # random.shuffle(doc_ids)
        # for doc_id in doc_ids[:10]:
        #     document = Document.objects.get(pk = doc_id)
        #     for section in document.section_set.all():
        #       section.validate = True
        #       section.save()


    def annotate_golden_documents(self, document_set):
        user, created = User.objects.get_or_create(username="goldenmaster")
        if created:
            user.set_password('')
            user.save()

        # Clean out all the old annotations just b/c we don't know what they were off on / need to be changed
        documents = Document.objects.filter(source = document_set).all()
        for doc in documents:
            views = View.objects.filter(section__document = doc, user = user)
            for view in views:
                Annotation.objects.filter(view = view).delete()

        with open('assets/datasets/'+ document_set +'_annos.txt','rU') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader, None)  # skip the headers
            for doc_id, doc_field, ann_type, text, start, stop in reader:
                try:
                  doc = Document.objects.get(document_id = doc_id)

                  for section in doc.section_set.all():
                      if section.kind == doc_field[0]:
                          view, created = View.objects.get_or_create(section = section, user = user)

                          ann, created = Annotation.objects.get_or_create(view = view, text = text, start = start, type = ann_type)
                          ann.kind = "e"
                          ann.user_agent = "goldenmaster"
                          ann.save()

                except Foo.DoesNotExist:
                    doc = None


        # Now go back over and confirm they match
        # gm_anns = Annotation.objects.filter(view__user = user).all()
        # for annotation in gm_anns:
        #   text = annotation.view.section.text
        #   print annotation.text, "::", text[annotation.start:]
        #   print "\n - - - - - - \n"




