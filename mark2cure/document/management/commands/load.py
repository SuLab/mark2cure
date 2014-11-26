from django.core.management.base import BaseCommand
from optparse import make_option

from mark2cure.document.models import Document, Section
from mark2cure.common.models import Task, DocumentQuestRelationship, SkillBadge

import random
import csv


class Command(BaseCommand):
    help = 'Load the initial documents into the database associate them with quests'

    option_list = BaseCommand.option_list + (
        make_option('--documents',
        action='store_true',
        dest='quest',
        default=False,
        help='Load the Documents into Mark2Cure and require 0 associations'),
    )

    option_list = BaseCommand.option_list + (
        make_option('--assign',
        action='store_true',
        dest='quest',
        default=False,
        help='Bin the Documents into assocaited Quests'),
    )

    option_list = BaseCommand.option_list + (
        make_option('--gm_anns',
        action='store_true',
        dest='quest',
        default=False,
        help='Load the Documents into Mark2Cure and require 0 associations'),
    )


    def handle(self, *args, **options):
        dataset = 'NCBI_corpus_testing_cleaned'
        print args
        print options

        if optionsp['gm_anns']:
            user, created = User.objects.get_or_create(username='Doc_G-man')
            if created:
                gold_profile = user.userprofile
                gold_profile.quote = "Be the change you wish to see in the world."
                gold_profile.save()

                user.set_password('')
                user.save()

                for index, level in enumerate(SkillBadge.levels):
                    badges.possibly_award_badge("skill_awarded", user=user, level=index+1)

                for task in Task.objects.all():
                    UserQuestRelationship.objects.create(task=task, user=user, completed=True)





            # Clean out all the old annotations just b/c we don't know what they were off on / need to be changed
            documents = Document.objects.filter(source=document_set).all()
            for doc in documents:
                views = View.objects.filter(section__document=doc, user=user)
                for view in views:
                    Annotation.objects.filter(view=view).delete()

            with open('assets/datasets/{0}_annos.txt'.format(document_set), 'rU') as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader, None)  # skip the headers
                for doc_id, doc_field, ann_type, text, start, stop in reader:
                    try:
                        doc = Document.objects.get(document_id=doc_id)

                        for section in doc.section_set.all():
                            if section.kind == doc_field[0]:
                                view, created = View.objects.get_or_create(section=section, user=user)
                                ann, created = Annotation.objects.get_or_create(view=view, text=text, start=start, type=ann_type)
                                ann.kind = 'e'
                                ann.user_agent = 'goldenmaster'
                                ann.save()

                    except DoesNotExist:
                        doc = None





        if options['assign']:
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
                for doc in task.documents.all():
                    dqr = DocumentQuestRelationship.objects.get(document=doc, task=task)
                    dqr.delete()

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

        if options['']:
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

