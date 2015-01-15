from django.core.management.base import BaseCommand
from optparse import make_option
from django.contrib.auth.models import User
from django.conf import settings

from mark2cure.userprofile.models import UserProfile
from mark2cure.document.models import Document, Section, View, Annotation
from mark2cure.common.models import Task, DocumentQuestRelationship, UserQuestRelationship, SkillBadge

from brabeion import badges

import random
import csv


class Command(BaseCommand):
    help = 'Utils for populating the system with Documents & Annotations'

    option_list = BaseCommand.option_list + (
        make_option('--documents',
            action='store_true',
            dest='documents',
            default=False,
            help='Load the Documents into Mark2Cure and require 0 associations'),

        make_option('--assign',
            action='store_true',
            dest='assign',
            default=False,
            help='Bin the Documents into associated Quests'),

        make_option('--gm_anns',
            action='store_true',
            dest='gm_anns',
            default=False,
            help='Load the Documents into Mark2Cure and require 0 associations'),
    )

    def handle(self, *args, **options):
        datasets = ['NCBI_corpus_testing', 'NCBI_corpus_training', 'NCBI_corpus_development']
        gm_documents_ids = [3464560, 7759075, 8198128, 3591825]

        doc_set = ['8589723', '8621452', '8644702', '8651278', '8661102',
                '8673131', '8675707', '8700509', '8751855', '8758207',
                '8786135', '8808605', '8824873', '8828602', '8843193',
                '8871666', '8929413', '8931709', '8944024', '8968716',
                '9012409', '9028321', '9056547', '9069115', '9090524',
                '9144439', '9174057', '8636252', '8689689', '8898652']
        mixed_gm = [8636252, 8689689, 8898652]

        '''
            Import all the annotations for our
            perceived "Expert" users. These annotations
            come from the GM Anns
        '''
        if options['gm_anns']:
            user, created = User.objects.get_or_create(username='Doc_G-man')
            if created:
                UserProfile.objects.create(user=user)

                gold_profile = user.userprofile
                gold_profile.quote = "Be the change you wish to see in the world."
                gold_profile.save()

                user.set_password('')
                user.save()

                # Assign the GM User the hightest skill possible
                for index, level in enumerate(SkillBadge.levels):
                    badges.possibly_award_badge("skill_awarded", user=user, level=index)

            # Clean out all the old annotations just b/c we don't know what they were off on / need to be changed
            for doc in Document.objects.all():
                views = View.objects.filter(section__document=doc, user=user)
                for view in views:
                    Annotation.objects.filter(view=view).delete()

            gm_documents_ids.extend(mixed_gm)
            for dataset in datasets:
                with open('assets/datasets/{dataset}_annos.txt'.format(dataset=dataset), 'r') as f:
                    reader = csv.reader(f, delimiter='\t')
                    next(reader, None)  # skip the headers
                    for doc_id, doc_field, ann_type, text, start, stop in reader:
                        if int(doc_id) in gm_documents_ids:

                            doc = Document.objects.get(document_id=doc_id)
                            print doc
                            dqr = DocumentQuestRelationship.objects.filter(document=doc).first()
                            print dqr.pk, dqr

                            # If the GM Anns are for a document that isn't included in any quests,
                            # don't bother
                            if dqr:
                                # Be (though uncomplete) associted with each first quest to link Views
                                gm_quest_rel, gm_quest_rel_created = UserQuestRelationship.objects.get_or_create(
                                    task=dqr.task,
                                    user=user)

                                for section in doc.section_set.all():
                                    # Make sure the annotion is for the title or abstract (our supported section types)
                                    if section.kind == doc_field[0]:
                                        view, created = View.objects.get_or_create(
                                            section=section,
                                            user=user,
                                            completed=True)

                                        gm_quest_rel.views.add(view)
                                        print gm_quest_rel.views.count()
                                        Annotation.objects.create(
                                            view=view,
                                            text=text,
                                            start=start,
                                            type=ann_type,
                                            kind='e')

        '''
            Randomly assign the loaded documents into our bins
            ensure the 4 training always go into Quest 1
        '''
        if options['assign']:
            gm_documents = Document.objects.filter(source='NCBI_corpus_training', document_id__in=gm_documents_ids).all()
            task_counter = 1

            # GM hidden
            task, task_created = Task.objects.get_or_create(
                name=str(task_counter),
                completions=None,
                requires_qualification=3,
                provides_qualification=4,
                points=5000,
                experiment=settings.EXPERIMENT)

            for doc in task.documents.all():
                dqr = DocumentQuestRelationship.objects.get(document=doc, task=task)
                dqr.delete()

            for gold_document in gm_documents:
                DocumentQuestRelationship.objects.create(task=task, document=gold_document)

            # Insert the other Documents
            document_set = list(
                Document.objects.filter(document_id__in=doc_set)
                                .exclude(document_id__in=gm_documents_ids)
                                .values_list('id', flat=True)
            )

            smallest_bin = 5
            largest_bin = 5
            completions = 3
            random.shuffle(document_set)

            while len(document_set) > smallest_bin:
                task_counter += 1
                task, task_created = Task.objects.get_or_create(
                    name=str(task_counter),
                    completions=completions,
                    requires_qualification=4,
                    provides_qualification=4,
                    points=5000,
                    experiment=settings.EXPERIMENT)

                quest_size = int(random.uniform(smallest_bin, largest_bin))
                sel = document_set[0:quest_size]

                for i in sel:
                    document_set.remove(i)

                for doc in task.documents.all():
                    dqr = DocumentQuestRelationship.objects.get(task=task, document=doc),
                    dqr.delete()

                for doc_idx in sel:
                    document = Document.objects.get(pk=doc_idx)
                    DocumentQuestRelationship.objects.create(task=task, document=document)

            else:
                if len(document_set) > 0:
                    task_counter += 1
                    task, task_created = Task.objects.get_or_create(
                        name=str(task_counter),
                        completions=completions,
                        requires_qualification=4,
                        provides_qualification=4,
                        points=5000,
                        experiment=settings.EXPERIMENT)

                    print "Adding", len(document_set), "to Quest:", task_counter, 'remaining', len(document_set)

                    for doc in task.documents.all():
                        dqr = DocumentQuestRelationship.objects.get(document=doc, task=task)
                        dqr.delete()

                    for doc_idx in document_set:
                        document = Document.objects.get(pk=doc_idx)
                        DocumentQuestRelationship.objects.create(task=task, document=document)

        '''
            Import the set of documents into Mark2Cure
            without performing any quest binning
        '''
        if options['documents']:
            for dataset in datasets:
                with open('assets/datasets/{dataset}_cleaned.txt'.format(dataset=dataset), 'r') as f:
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

