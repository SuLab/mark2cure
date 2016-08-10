from django.core.urlresolvers import reverse
from django.test import TestCase

from ..task.models import Task, UserQuestRelationship
from ..common.models import Group
from .tasks import get_pubmed_document
from .models import Document, Section, Pubtator, Annotation
from ..test_base.test_base import TestBase

from ..common.bioc import BioCReader
from Bio import Entrez
import datetime
import json


class DocumentImportProcessing(TestCase):

    def setUp(self):
        # random PMID from NCBI_corpus_testing
        self.pmid = 9467011
        get_pubmed_document(self.pmid)
        self.doc = Document.objects.get(document_id=self.pmid)

    def test_document_init(self):
        # Test the document loaded to DB
        self.assertEqual(Document.objects.count(), 1)
        self.assertEqual(self.doc.document_id, self.pmid)

        # Test the sections were properly included
        sections = Section.objects.filter(document=self.doc)
        self.assertEqual(sections.count(), 2)

        # Test pubtator calls were made
        # (TODO) Actually monitor the session-id of these requests
        pubtators = Pubtator.objects.filter(document=self.doc)
        self.assertEqual(pubtators.count(), 3)

        # Test document specific methods
        self.assertEqual(self.doc.available_sections().count(), 2)
        self.assertEqual(self.doc.count_available_sections(), 2)

        # Actually determine if a doc would be repadded
        self.assertEqual(self.doc.update_padding(), False)


class PubtatorImportProcessing(TestCase):

    def setUp(self):
        date = datetime.datetime.today() - datetime.timedelta(days=5)
        h = Entrez.esearch(db='pubmed', retmax=10, term='("{date}"[Date - Publication] : "3000"[Date - Publication])'.format(date=date.strftime('%Y/%m/%M')))
        result = Entrez.read(h)
        for pmid in result.get('IdList'):
            print pmid

    def test_document_init(self):
        pass


class DocumentAPIMethods(TestCase):
    pass


class DocumentAPIViews(TestCase):
    fixtures = ['tests_document.json']

    @classmethod
    def setUp(cls):
        cls.doc = Document.objects.first()

    def test_document_init(self):
        self.assertEqual(Document.objects.count(), 1)
        self.assertEqual(Section.objects.count(), 2)
        self.assertEqual(Pubtator.objects.count(), 3)

    def test_document_as_bioc(self):
        response = self.client.get('/document/{pmid}.json'.format(pmid=self.doc.document_id))
        json_string = response.content
        self.assertNotEqual(json_string, '', msg='API returned empty response for document BioC Representation.')

        data = json.loads(json_string)

        # Make sure it's the same document in BioC as DB
        self.assertEqual(int(data.get('collection').get('document').get('id')),
                         self.doc.document_id)
        self.assertEqual(len(data.get('collection').get('document').get('passage')), 2)
        self.assertEqual(data.get('collection').get('document').get('passage')[0].get('text'),
                         self.doc.section_set.first().text)
        self.assertEqual(data.get('collection').get('document').get('passage')[1].get('text'),
                         self.doc.section_set.last().text)

        # Make sure it doesn't contain any annotations
        self.assertEqual(data.get('collection').get('document').get('passage')[0].get('annotation'), [None])
        self.assertEqual(data.get('collection').get('document').get('passage')[1].get('annotation'), [None])

        # We already validated everything in JSON b/c it's easier. Let's just
        # make sure the XML document passes too without specific checks
        response = self.client.get('/document/{pmid}.xml'.format(pmid=self.doc.document_id))
        r = BioCReader(source=response.content)
        r.read()
        self.assertEqual(len(r.collection.documents), 1)
        self.assertEqual(int(r.collection.documents[0].id), self.doc.document_id)
        self.assertEqual(len(r.collection.documents[0].passages), 2)
        self.assertEqual(len(r.collection.documents[0].passages[0].annotations), 0)
        self.assertEqual(len(r.collection.documents[0].passages[1].annotations), 0)

    def test_document_as_bioc_with_pubtator(self):
        # Unused for now
        """
        pub_query_set = Pubtator.objects.filter(
                document=self.doc,
                session_id='',
                content__isnull=False)
        """
        response = self.client.get('/document/pubtator/{pmid}.json'.format(pmid=self.doc.document_id))
        json_string = response.content
        self.assertNotEqual(json_string, '', msg='API returned empty response for document BioC Pubtator Representation.')

        data = json.loads(json_string)

        # Make sure it's the same document in BioC as DB
        self.assertEqual(int(data.get('collection').get('document').get('id')), self.doc.document_id)
        self.assertEqual(len(data.get('collection').get('document').get('passage')), 2)
        self.assertEqual(data.get('collection').get('document').get('passage')[0].get('text'), self.doc.section_set.first().text)
        self.assertEqual(data.get('collection').get('document').get('passage')[1].get('text'), self.doc.section_set.last().text)

        # Make sure it contains any annotations
        self.assertNotEqual(len(data.get('collection').get('document').get('passage')[0].get('annotation')), 0)
        self.assertNotEqual(len(data.get('collection').get('document').get('passage')[1].get('annotation')), 0)

        # We already validated everything in JSON b/c it's easier. Let's just
        # make sure the XML document passes too without specific checks
        response = self.client.get('/document/pubtator/{pmid}.xml'.format(pmid=self.doc.document_id))
        r = BioCReader(source=response.content)
        r.read()
        self.assertEqual(len(r.collection.documents), 1)
        self.assertEqual(int(r.collection.documents[0].id), self.doc.document_id)
        self.assertEqual(len(r.collection.documents[0].passages), 2)
        self.assertNotEqual(len(r.collection.documents[0].passages[0].annotations), 0)
        self.assertNotEqual(len(r.collection.documents[0].passages[1].annotations), 0)


class DocumentSubmissionsAPIViews(TestCase, TestBase):
    fixtures = ['tests_documents.json', 'tests_common.json']

    @classmethod
    def setUp(cls):

        cls.group = Group.objects.first()
        cls.task = Task.objects.first()

        cls.user_names = ['player', 'opponent']

        cls.users = {}
        cls.user_annotation_list = []

    def test_document_as_bioc_for_pairing(self):
        self.create_new_user_accounts(self.user_names)

        # Ensure the player views the Q but can't match b/c no Anns exist
        self.client.login(username='player', password='password')

        # Ensure the User info is showing up in the header
        response = self.client.get('/dashboard/')
        self.assertInHTML('<p>Level: Expert</p>', response.content)

        # Ensure no User >> Quest views until after viewed once
        self.assertEqual(UserQuestRelationship.objects.count(), 0)
        response = self.client.get(reverse('task-entity-recognition:quest-home',
                                   kwargs={'quest_pk': self.task.pk}),
                                   follow=True)
        doc = self.get_document_from_response(self.user_names[0])
        self.assertEqual(UserQuestRelationship.objects.count(), 1)

        # Ensure this returns a 500 for the player b/c there are no submissions yet
        response = self.client.get(reverse('task-entity-recognition:results-bioc',
                                           kwargs={'task_pk': self.task.pk,
                                           'doc_pk': doc.pk,
                                           'format_type': 'xml'}))
        self.assertEqual(response.status_code, 500)
        # (TODO) this test broke, response.content = "points awarded" now?
        self.assertEqual(response.content, 'no_points_awarded')
        self.client.logout()
        # Submit bogus Annotations as opponent to try match again for player
        self.client.login(username='opponent', password='password')
        self.assertEqual(Annotation.objects.count(), 0)
        response = self.client.get(reverse('task-entity-recognition:quest-home',
                                   kwargs={'quest_pk': self.task.pk}),
                                   follow=True)

        # load fake annotations from method in TestBase, where this checks for
        # correct response status.
        self.load_fake_annotations(self.user_names, doc)

        # There should be many annotations after loading fake annotations
        self.assertGreater(Annotation.objects.count(), 3)

        # Then submit the document for the Quest
        response = self.client.post(reverse('task-entity-recognition:doc-quest-submit',
                                    kwargs={'quest_pk': self.task.pk,
                                    'document_pk': doc.pk}), follow=True)
        self.client.logout()

        # Try again as the player to see if comparison uses opponents
        self.client.login(username='player', password='password')

        # Submit this Document without contributing any Annotations
        response = self.client.post(reverse('task-entity-recognition:doc-quest-submit',
                                    kwargs={'quest_pk': self.task.pk,
                                    'document_pk': doc.pk}), follow=True)
        # Fetch the BioC Document again
        response = self.client.get(reverse('task:results-bioc',
                                           kwargs={'task_pk': self.task.pk,
                                           'doc_pk': doc.pk,
                                           'format_type': 'xml'}))
        self.assertEqual(response.status_code, 200)
        r = BioCReader(source=response.content)
        r.read()

        # Make sure the BioC document has the opponent's info
        self.assertEqual(len(r.collection.documents), 1)
        self.assertEqual(int(r.collection.documents[0].id), doc.document_id)
        self.assertEqual(len(r.collection.documents[0].passages), 2)

        # 1st passage with 2 players x 30 maximum annotations, gives range:
        self.assertTrue(0 <= len(r.collection.documents[0].passages[0].annotations) <= 60)

        # 2nd passage with 2 players with 30 max annotations/player, gives range:
        self.assertTrue(0 <= len(r.collection.documents[0].passages[1].annotations) <= 60)

        self.assertEqual(r.collection.documents[0].passages[1].annotations[0].infons['user_name'], 'opponent')

        # Annotation types are randomized now
        annotation_type = int(r.collection.documents[0].passages[1].annotations[0].infons['type'])
        self.assertTrue(annotation_type == 0 or annotation_type == 1 or annotation_type == 2)

        # Generates random text, so cannot test exact annotations, so check for str entry
        self.assertEqual(type(r.collection.documents[0].passages[1].annotations[0].text), type('string'))
        self.assertEqual(int(r.collection.infons['points']), 0)
        self.assertEqual(r.collection.infons['partner'], 'opponent')
        self.client.logout()

    def test_document_as_bioc_with_m2c(self):
        # Submit Annotations (As User 1) so they show up when inspecting the M2C submissions
        self.create_new_user_accounts(self.user_names)
        self.assertEqual(Annotation.objects.count(), 0)
        self.client.login(username='player', password='password')
        response = self.client.get(reverse('task-entity-recognition:quest-home',
                                   kwargs={'quest_pk': self.task.pk}),
                                   follow=True)

        doc = response.context['document']
        abstract = doc.available_sections().last()

        # Annotation submit URL
        url = reverse('task-entity-recognition:create', kwargs={'task_pk': int(self.task.pk),
                      'section_pk': int(abstract.pk)})
        self.assertEqual(self.client.post(url, {'type': 0,
                                          'text': 'text annotation 0',
                                          'start': 0}).status_code, 200)
        self.assertEqual(self.client.post(url, {'type': 1,
                                          'text': 'text annotation 1',
                                          'start': 10}).status_code, 200)
        self.assertEqual(self.client.post(url, {'type': 2,
                                          'text': 'text annotation 2',
                                          'start': 20}).status_code, 200)
        self.assertEqual(Annotation.objects.count(), 3)

        # Then submit the document for the Quest
        response = self.client.post(reverse('task-entity-recognition:doc-quest-submit',
                                    kwargs={'quest_pk': self.task.pk,
                                    'document_pk': doc.pk}),
                                    follow=True)
        self.client.logout()

        # Submit Annotations (As User 1) so they show up when inspecting the M2C submissions
        self.assertEqual(Annotation.objects.count(), 3)
        self.client.login(username='opponent', password='password')
        response = self.client.get(reverse('task-entity-recognition:quest-home',
                                   kwargs={'quest_pk': self.task.pk}),
                                   follow=True)

        # Annotation submit URL
        url = reverse('task-entity-recognition:create', kwargs={'task_pk': self.task.pk,
                      'section_pk': abstract.pk})
        self.assertEqual(self.client.post(url, {'type': 0,
                                          'text': 'text annotation 3',
                                          'start': 30}).status_code, 200)
        self.assertEqual(self.client.post(url, {'type': 1,
                                          'text': 'text annotation 4',
                                          'start': 40}).status_code, 200)
        self.assertEqual(self.client.post(url, {'type': 2,
                                          'text': 'text annotation 5',
                                          'start': 50}).status_code, 200)
        self.assertEqual(Annotation.objects.count(), 6)

        # Then submit the document for the Quest
        response = self.client.post(reverse('task-entity-recognition:doc-quest-submit',
                                    kwargs={'quest_pk': self.task.pk,
                                    'document_pk': doc.pk}),
                                    follow=True)
        self.client.logout()

