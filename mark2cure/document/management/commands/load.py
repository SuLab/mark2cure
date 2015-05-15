from django.core.management.base import BaseCommand
from optparse import make_option
from django.contrib.auth.models import User
from django.conf import settings

from mark2cure.userprofile.models import UserProfile
from mark2cure.document.models import Document, Pubtator, Section, View, Annotation

from mark2cure.common.models import Task, Group, DocumentQuestRelationship, UserQuestRelationship, SkillBadge
from mark2cure.document.tasks import get_pubmed_document, get_pubtator_response

from brabeion import badges

import requests
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

        make_option('--pubtator',
            action='store_true',
            dest='pubtator',
            default=False,
            help='Fetch any lingering pubtator content'),

        make_option('--assign',
            action='store_true',
            dest='assign',
            default=False,
            help='Bin the Documents into associated Quests'),

        make_option('--training',
            action='store_true',
            dest='training',
            default=False,
            help='Load the Training Documents into Mark2Cure'),

        make_option('--gm_anns',
            action='store_true',
            dest='gm_anns',
            default=False,
            help='Load the Documents into Mark2Cure and require 0 associations'),
    )

    def handle(self, *args, **options):
        datasets = ['NCBI_corpus_testing', 'NCBI_corpus_training', 'NCBI_corpus_development']
        gm_documents_ids = [3464560, 7759075, 8198128, 3591825]

        # doc_set = ['8589723', '8621452', '8644702', '8651278', '8661102',
        #        '8673131', '8675707', '8700509', '8751855', '8758207',
        #        '8786135', '8808605', '8824873', '8828602', '8843193',
        #        '8871666', '8929413', '8931709', '8944024', '8968716',
        #        '9012409', '9028321', '9056547', '9069115', '9090524',
        #        '9144439', '9174057', '8636252', '8689689', '8898652']
        # mixed_gm = [8636252, 8689689, 8898652]

        doc_set = Document.objects.filter(source='NCBI_corpus_training').values_list('document_id', flat=True)
        # 589 * .1 = 58.9
        mixed_gm = list(doc_set)
        random.shuffle(mixed_gm)
        mixed_gm = mixed_gm[:59]

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
                            dqr = DocumentQuestRelationship.objects.filter(document=doc, task__experiment=settings.EXPERIMENT).first()
                            print "Document ID:", doc.pk, "Quest Relationship:", dqr.pk
                            return

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
                                        Annotation.objects.create(
                                            view=view,
                                            text=text,
                                            start=start,
                                            type=ann_type,
                                            kind='e')

        '''
            1) Create a Group to Organize Tasks
            2) Create tasks w/ 5 documents each that all have valid pubtator
                content from 622 set
            3) Organize these tasks into Group to know wh
        '''
        if options['assign']:
            '''
            # GM hidden
            task, task_created = Task.objects.get_or_create(
                name=str(task_counter),
                completions=None,
                requires_qualification=,
                provides_qualification=4,
                points=5000,
                experiment=settings.EXPERIMENT)

            for doc in task.documents.all():
                dqr = DocumentQuestRelationship.objects.get(document=doc, task=task)
                dqr.delete()

            for gold_document in gm_documents:
                DocumentQuestRelationship.objects.create(task=task, document=gold_document)
            '''

            documents = Document.objects.filter(source='pubmed')
            task_counter = Task.objects.filter(kind='q').count()

            # Insert the other Documents
            document_set = list(documents.values_list('id', flat=True))

            smallest_bin = 5
            largest_bin = 5
            completions = 15
            random.shuffle(document_set)

            group, group_v = Group.objects.get_or_create(name='622', stub='622', enabled=True)
            #task.clear_documents()
            last_task = group.task_set.last()

            while len(document_set) > smallest_bin:

                quest_size = int(random.uniform(smallest_bin, largest_bin))
                # If there was an existing Task with less than the
                # desired number of documents
                if last_task and last_task.documents.count() < quest_size:

                    # Shuffle & Remove the document_pk for use and from being selected again
                    random.shuffle(document_set)
                    doc_pk = document_set[0]
                    document_set.remove(doc_pk)

                    document = Document.objects.get(pk=doc_pk)
                    print 'Add Document', len(document_set), document.valid_pubtator(), last_task
                    if document.valid_pubtator():
                        DocumentQuestRelationship.objects.create(task=last_task, document=document)

                else:
                    print '> Add New Task'
                    if last_task:
                        idx = last_task.pk
                    else:
                        idx = Task.objects.last().pk

                    last_task, task_created = Task.objects.get_or_create(
                        name=str(idx+1),
                        completions=completions,
                        requires_qualification=7,
                        provides_qualification=7,
                        points=5000,
                        group=group)


        if options['training']:
            group, group_c = Group.objects.get_or_create(name='Training', stub='training', enabled=False)

            # Disease marking
            for pmid in [25514328, 21431621, 10364520]:
                get_pubmed_document(pmid, 'training-disease')
                task, task_created = Task.objects.get_or_create(name='T1', completions=None, requires_qualification=7, provides_qualification=7, points=5000, group=group)
                document = Document.objects.get(document_id=pmid, source='training-disease', include_pubtator=False)
                DocumentQuestRelationship.objects.get_or_create(task=task, document=document)

            # Genes marking
            for pmid in [23542699, 25274141, 23437350]:
                get_pubmed_document(pmid, 'training-genes')
                task, task_created = Task.objects.get_or_create(name='T2', completions=None, requires_qualification=7, provides_qualification=7, points=5000, group=group)
                document = Document.objects.get(document_id=pmid, source='training-genes', include_pubtator=False)
                DocumentQuestRelationship.objects.get_or_create(task=task, document=document)

            # Treatment marking
            for pmid in [25467138, 25671401, 25779362]:
                get_pubmed_document(pmid, 'training-treatment')
                task, task_created = Task.objects.get_or_create(name='T3', completions=None, requires_qualification=7, provides_qualification=7, points=5000, group=group)
                document = Document.objects.get(document_id=pmid, source='training-treatment', include_pubtator=False)
                DocumentQuestRelationship.objects.get_or_create(task=task, document=document)

            # All concepts
            for pmid in [24883236, 25732996, 18797263, 25663566]:
                get_pubmed_document(pmid, 'training-all')
                task, task_created = Task.objects.get_or_create(name='T4', completions=None, requires_qualification=7, provides_qualification=7, points=5000, group=group)
                document = Document.objects.get(document_id=pmid, source='training-all', include_pubtator=False)
                DocumentQuestRelationship.objects.get_or_create(task=task, document=document)

            # Extra practice
            for pmid in [25927578, 25931357]:
                get_pubmed_document(pmid, 'training-extra')
                task, task_created = Task.objects.get_or_create(name='T5', completions=None, requires_qualification=7, provides_qualification=7, points=5000, group=group)
                document = Document.objects.get(document_id=pmid, source='training-extra', include_pubtator=False)
                DocumentQuestRelationship.objects.get_or_create(task=task, document=document)


        '''
            Import the set of documents into Mark2Cure
            without performing any quest binning
        '''
        if options['documents']:
            # https://s3.amazonaws.com/uploads.hipchat.com/25885/154162/RD9pzKYNggkhj6O/combined.txt
            res = requests.get('https://s3.amazonaws.com/uploads.hipchat.com/25885/154162/fVP7w1pKOOQCOJa/congen_dis_glyco_pmids.txt')
            ids = res.text.split('\n')
            for pmid in ids:
                get_pubmed_document(pmid)

        '''
            Check for any unfetched Pubtator objects
        '''
        if options['pubtator']:
            '''
            from django.db.models.signals import post_save
            for pubtator in Pubtator.objects.filter(content__isnull=True).all():
                post_save.send(Pubtator, instance=pubtator, created=True)
            '''

            '''
            # Gracefully check for updates to session id
            for pubtator in Pubtator.objects.filter(content__isnull=True).all():
                # Make response to post job to pubtator
                payload = {'content-type': 'text/xml'}
                writer = pubtator.document.as_writer()
                data = str(writer)

                get_pubtator_response.apply_async(
                    args=[pubtator.pk, data, payload, 0],
                )
            '''

