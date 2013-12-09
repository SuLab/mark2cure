from mark2cure.document.models import Document, Section
from mark2cure.document.utils import strip_tags
from django.core.exceptions import ObjectDoesNotExist

import csv

path = "NCBI_corpus_development"
with open('assets/NCBI_corpus/'+ path +'.txt','r') as f:
    reader = csv.reader(f, delimiter='\t')
    for num, title, text in reader:
        title = title.replace('<category=', '<category type=')
        plain_title = str(strip_tags(title, ['category']))

        text = text.replace('<category=', '<category type=')
        plain_text = str(strip_tags(text, ['category']))

        print plain_title

        try:
            doc = Document.objects.get(document_id = num)
        except ObjectDoesNotExist:
            doc = Document()

            doc.document_id = num
            doc.title = plain_title
            doc.source = path
            doc.save()

            sec = Section(kind = "t")
            sec.text = plain_title
            sec.document = doc
            sec.save()

            sec = Section(kind = "a")
            sec.text = plain_text
            sec.document = doc
            sec.save()


