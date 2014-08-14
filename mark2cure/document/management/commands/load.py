from django.core.management.base import BaseCommand

from mark2cure.document.models import Document, Section

import csv


class Command(BaseCommand):

    def handle(self, *args, **options):
        command = args[0]
        self.stdout.write('-- Loading Documents: ({0}) --'.format(command))

        with open('assets/datasets/{0}-padded.txt'.format(command), 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for num, title, text in reader:
                print title

                doc, doc_c = Document.objects.get_or_create(document_id=num)
                doc.title = title
                doc.source = command
                doc.save()

                sec, sec_c = Section.objects.get_or_create(kind='t', document=doc)
                sec.text = title
                sec.save()

                sec, sec_c = Section.objects.get_or_create(kind='a', document=doc)
                sec.text = text
                sec.save()




