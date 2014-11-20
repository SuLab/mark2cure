from django.core.management.base import BaseCommand
from optparse import make_option

from mark2cure.document.models import Document, Section
from mark2cure.common.models import Task, DocumentQuestRelationship

import random
import csv


class Command(BaseCommand):
    help = 'Load the initial documents into the database associate them with quests'

    option_list = BaseCommand.option_list + (
        make_option('--quest',
        action='store_true',
        dest='quest',
        default=False,
        help='Bin the Documents into assocaited Quests'),
    )

    def handle(self, *args, **options):
        dataset = 'NCBI_corpus_testing_cleaned'
        print args
        print options

        if options['quest']:
            documents = list(Document.objects.filter(source=dataset).values_list('id', flat=True))
            smallest_bin = 3
            largest_bin = 10
            random.shuffle(documents)
            quest_id = 3

            while len(documents) > smallest_bin:
                quest_size = int(random.uniform(smallest_bin, largest_bin))
                sel = documents[0:quest_size]

                for i in sel: documents.remove(i)
                print "Adding", len(sel), "to Quest:", quest_id
                task = Task.objects.get(pk=quest_id)
                for doc in task.documents.all(): task.documents.remove(doc)

                for i in sel:
                    document = Document.objects.get(pk=i)
                    DocumentQuestRelationship.objects.create(task=task, document=document)

                quest_id += 1
            else:
                if len(documents) > 0:
                    print "Adding", len(documents), "to Quest:", quest_id
                    task = Task.objects.get(pk=quest_id)
                    for doc in task.documents.all(): task.documents.remove(doc)

                    for i in sel:
                        document = Document.objects.get(pk=i)
                        DocumentQuestRelationship.objects.create(task=task, document=document)

        else:
            with open('assets/datasets/{dataset}.txt'.format(dataset=dataset), 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for num, title, text in reader:
                    print title
                    doc, doc_c = Document.objects.get_or_create(document_id=num)
                    doc.title = title
                    doc.source = dataset
                    doc.save()

                    sec, sec_c = Section.objects.get_or_create(kind='t', document=doc)
                    sec.text = title
                    sec.save()

                    sec, sec_c = Section.objects.get_or_create(kind='a', document=doc)
                    sec.text = text
                    sec.save()

