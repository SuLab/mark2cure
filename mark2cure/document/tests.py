from django.test import TestCase
from rest_framework.test import APIRequestFactory
from django.core.urlresolvers import reverse

from .tasks import get_pubmed_document
from .models import Document, Section, Pubtator, Annotation, View
from ..common.models import Group, Task, UserQuestRelationship
from django.contrib.auth.models import User
from brabeion import badges

from mark2cure.common.bioc import BioCReader
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
        self.assertEqual(int(data.get('collection').get('document').get('id')), self.doc.document_id)
        self.assertEqual(len(data.get('collection').get('document').get('passage')), 2)
        self.assertEqual(data.get('collection').get('document').get('passage')[0].get('text'), self.doc.section_set.first().text)
        self.assertEqual(data.get('collection').get('document').get('passage')[1].get('text'), self.doc.section_set.last().text)

        # Make sure it doesn't contain any annotations
        self.assertEqual(data.get('collection').get('document').get('passage')[0].get('annotation'), None)
        self.assertEqual(data.get('collection').get('document').get('passage')[1].get('annotation'), None)

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
        pub_query_set = Pubtator.objects.filter(
                document=self.doc,
                session_id='',
                content__isnull=False)

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


class DocumentSubmissionsAPIViews(TestCase):
    fixtures = ['tests_documents.json', 'tests_common.json']

    @classmethod
    def setUp(cls):
        cls.group = Group.objects.first()
        cls.task = Task.objects.first()

        cls.player = User.objects.create_user('player', password='password')
        badges.possibly_award_badge("skill_awarded", user=cls.player, level=7, force=True)

        cls.opponent = User.objects.create_user('opponent', password='password')
        badges.possibly_award_badge("skill_awarded", user=cls.opponent, level=7, force=True)

    def test_document_as_bioc_for_pairing(self):
        # Ensure the player views the Q but can't match b/c no Anns exist
        self.client.login(username='player', password='password')

        # Ensure the User info is showing up in the header
        response = self.client.get('/dashboard/')
        self.assertInHTML('<p>Level: Expert</p>', response.content)

        # Ensure no User >> Quest views until after viewed once
        self.assertEqual(UserQuestRelationship.objects.count(), 0)
        response = self.client.get(reverse('common:quest-home', kwargs={'quest_pk': self.task.pk}), follow=True)
        doc = response.context['document']
        self.assertEqual(UserQuestRelationship.objects.count(), 1)

        # Ensure this returns a 500 for the player b/c there are no submissions yet
        response = self.client.get(reverse('document:results-bioc', kwargs={'task_pk': self.task.pk, 'doc_pk': doc.pk, 'format_type': 'xml'}))
        self.assertEqual(response.status_code, 500)
        self.assertEqual(response.content, 'no_points_awarded')
        self.client.logout()

        #
        # Submit bogus Annotations as opponent to try match again for player
        #
        self.client.login(username='opponent', password='password')
        self.assertEqual(Annotation.objects.count(), 0)
        response = self.client.get(reverse('common:quest-home', kwargs={'quest_pk': self.task.pk}), follow=True)

        # Annotation submit URL
        abstract = doc.available_sections().last()
        url = reverse('document:create', kwargs={'task_pk': self.task.pk, 'section_pk': abstract.pk})
        self.client.post(url, {'type': 0, 'text': 'text annotation 0', 'start': 0})
        self.client.post(url, {'type': 1, 'text': 'text annotation 1', 'start': 10})
        self.client.post(url, {'type': 2, 'text': 'text annotation 2', 'start': 20})
        self.assertEqual(Annotation.objects.count(), 3)

        # Then submit the document for the Quest
        response = self.client.post(reverse('common:doc-quest-submit', kwargs={'quest_pk': self.task.pk, 'document_pk': doc.pk}), follow=True)
        self.client.logout()


        #
        # Try again as the player to see if comparison uses opponents
        #
        self.client.login(username='player', password='password')
        # Submit this Document without contributing any Annotations
        response = self.client.post(reverse('common:doc-quest-submit', kwargs={'quest_pk': self.task.pk, 'document_pk': doc.pk}), follow=True)

        # Fetch the BioC Document again
        response = self.client.get(reverse('document:results-bioc', kwargs={'task_pk': self.task.pk, 'doc_pk': doc.pk, 'format_type': 'xml'}))
        self.assertEqual(response.status_code, 200)
        r = BioCReader(source=response.content)
        r.read()

        self.assertEqual(len(r.collection.documents), 1)
        self.assertEqual(int(r.collection.documents[0].id), doc.document_id)
        self.assertEqual(len(r.collection.documents[0].passages), 2)
        self.assertEqual(len(r.collection.documents[0].passages[0].annotations), 0)
        self.assertEqual(len(r.collection.documents[0].passages[1].annotations), 3)
        self.client.logout()


'''
    def test_document_as_bioc_with_m2c(self):
        pass
'''

