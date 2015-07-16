from django.test import TestCase
from rest_framework.test import APIRequestFactory
from django.core.urlresolvers import reverse

from .tasks import get_pubmed_document
from .models import Document, Section, Pubtator, Annotation, View
from ..common.models import Group, Task, UserQuestRelationship
from ..test_base.test_base import TestBase
from django.contrib.auth.models import User
from brabeion import badges

from mark2cure.common.bioc import BioCReader
from Bio import Entrez, Medline
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
        # TODO max question.  prints ids no matter what so how to test it?  Pubtator isn't working now and it is printing IDs.  I am confused about this.
        for i in pubtators:
            print i.session_id

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
        pub_query_set = Pubtator.objects.filter(document=self.doc,
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


class DocumentSubmissionsAPIViews(TestCase, TestBase):
    """class DocumentSubmissionsAPIViews can now inherit methods from TestBase,
    which is a class located in mark2cure/test_base/test_base.py
    """
    # snapshot of a specific table.
    fixtures = ['tests_documents.json', 'tests_common.json']

    @classmethod
    def setUp(cls):
        cls.group = Group.objects.first()
        cls.task = Task.objects.first()

        cls.user_names = ['UserA', 'UserB']
        cls.users = {}
        cls.user_annotation_list = []

        for user_name in cls.user_names:
            cls.users[user_name] = User.objects.create_user(user_name, password='password')
            badges.possibly_award_badge("skill_awarded", user=cls.users[user_name], level=7, force=True)

    def test_document_as_bioc_for_pairing(self):
        # Ensure the player views the Q but can't match b/c no Anns exist TODO # Jennifer
        # just use the first user to test "viewing the quest"
        user_name = self.user_names[0]
        self.client.login(username=user_name, password='password')

        # Ensure the User info is showing up in the header
        response = self.client.get('/dashboard/')
        self.assertInHTML('<p>Level: Expert</p>', response.content)

        # Ensure no User >> Quest views until after viewed once
        self.assertEqual(UserQuestRelationship.objects.count(), 0)
        # this happens one time
        response = self.client.get(reverse('common:quest-home', kwargs={'quest_pk': self.task.pk}), follow=True)
        doc = response.context['document']
        # Ensure one view after Quest views
        self.assertEqual(UserQuestRelationship.objects.count(), 1)

        # Ensure this returns a 500 for the player b/c there are no submissions yet
        response = self.client.get(reverse('document:results-bioc', kwargs={'task_pk': self.task.pk, 'doc_pk': doc.pk, 'format_type': 'xml'}))
        # TODO: max question, does the user get 500 pts even if they don't do any submissions?  No disease,gene,drugs?
        self.assertEqual(response.status_code, 500)
        self.assertEqual(response.content, 'no_points_awarded')
        self.client.logout()

        # Submit bogus Annotations as opponent to try match again for player
        opponent = self.user_names[1]  # TODO use more opponents
        self.client.login(username=opponent, password='password')
        # print Annotation.objects.count()  # testing to see how Annotation.objects.count() looks
        self.assertEqual(Annotation.objects.count(), 0)  # Need to set annotations to 0 again
        self.client.logout()

        opponent = self.user_names[1]  # TODO use more opponents
        self.client.login(username=opponent, password='password')
        # load fake annotations using method in TestBase class
        self.load_fake_annotations()
        self.client.get(reverse('common:quest-home', kwargs={'quest_pk': self.task.pk}), follow=True)
        # Annotation submit URL TODO (this is already done in the load fake annotations, no need to repeat it)

        """
        abstract = doc.available_sections().last()
        url = reverse('document:create', kwargs={'task_pk': self.task.pk, 'section_pk': abstract.pk})
        self.assertEqual(self.client.post(url, {'type': 0, 'text': 'text annotation 0', 'start': 0}).status_code, 200)
        self.assertEqual(self.client.post(url, {'type': 1, 'text': 'text annotation 1', 'start': 10}).status_code, 200)
        self.assertEqual(self.client.post(url, {'type': 2, 'text': 'text annotation 2', 'start': 20}).status_code, 200)
        """
        # TODO there should be many annotations (range?  Check this value)
        # return the integer for the number of annotations
        self.assertGreater(Annotation.objects.count(), 3)

        # Then submit the document for the Quest
        response = self.client.post(reverse('common:doc-quest-submit', kwargs={'quest_pk': self.task.pk, 'document_pk': doc.pk}), follow=True)
        self.client.logout()

        # Try again as the player to see if comparison uses opponents
        user_name = self.user_names[0]  # Jennifer TODO
        self.client.login(username=user_name, password='password')
        # Submit this Document without contributing any Annotations
        response = self.client.post(reverse('common:doc-quest-submit', kwargs={'quest_pk': self.task.pk, 'document_pk': doc.pk}), follow=True)

        # Fetch the BioC Document again
        response = self.client.get(reverse('document:results-bioc', kwargs={'task_pk': self.task.pk, 'doc_pk': doc.pk, 'format_type': 'xml'}))
        self.assertEqual(response.status_code, 200) #  TODO

        r = BioCReader(source=response.content)
        r.read()

        # Make sure the BioC document has the opponent's info
        opponent = self.user_names[1]  # Jennifer TODO
        self.assertEqual(len(r.collection.documents), 1)
        self.assertEqual(int(r.collection.documents[0].id), doc.document_id)
        self.assertEqual(len(r.collection.documents[0].passages), 2)

        # assert that there are 0 to 30 annotations for passage 1 and 2
        self.assertTrue(0 <= len(r.collection.documents[0].passages[0].annotations) <= 30)
        self.assertTrue(0 <= len(r.collection.documents[0].passages[1].annotations) <= 30)

        # TODO, ask Max about infons (sending me the paper on bioC...  artifact of file format... add metadata to a response, use an infon to response, (score of partner, and lev))

        # self.assertEqual(r.collection.documents[0].passages[1].annotations[0].infons['user_name'], user_name)  # TODO fix this
        # self.assertEqual(int(r.collection.documents[0].passages[1].annotations[0].infons['type']), 0)  # TODO
        # self.assertEqual(r.collection.documents[0].passages[1].annotations[0].text, 'text annotation 0')  # TODO
        # you don't submit any annotations and your partner has some, you get zero points
        self.assertEqual(int(r.collection.infons['points']), 0)
        # test when there are multiple people annotating # TODO
        self.assertEqual(r.collection.infons[opponent], 'UserB') #
        self.client.logout()

    def test_document_as_bioc_with_m2c(self):
        # TODO Max: these very basic text annotations are okay because we can test the more sophisticated ones later
        # Submit Annotations (As User 1) so they show up when inspecting the M2C submissions
        self.assertEqual(Annotation.objects.count(), 0)
        self.client.login(username='player', password='password')
        response = self.client.get(reverse('common:quest-home', kwargs={'quest_pk': self.task.pk}), follow=True)

        doc = response.context['document']
        abstract = doc.available_sections().last()

        # Annotation submit URL
        url = reverse('document:create', kwargs={'task_pk': self.task.pk, 'section_pk': abstract.pk})
        self.assertEqual(self.client.post(url, {'type': 0, 'text': 'text annotation 0', 'start': 0}).status_code, 200)
        self.assertEqual(self.client.post(url, {'type': 1, 'text': 'text annotation 1', 'start': 10}).status_code, 200)
        self.assertEqual(self.client.post(url, {'type': 2, 'text': 'text annotation 2', 'start': 20}).status_code, 200)
        self.assertEqual(Annotation.objects.count(), 3)
        # Then submit the document for the Quest
        response = self.client.post(reverse('common:doc-quest-submit', kwargs={'quest_pk': self.task.pk, 'document_pk': doc.pk}), follow=True)
        self.client.logout()

        # Submit Annotations (As User 1) so they show up when inspecting the M2C submissions
        self.assertEqual(Annotation.objects.count(), 3)
        self.client.login(username='opponent', password='password')
        response = self.client.get(reverse('common:quest-home', kwargs={'quest_pk': self.task.pk}), follow=True)

        # Annotation submit URL
        url = reverse('document:create', kwargs={'task_pk': self.task.pk, 'section_pk': abstract.pk})
        self.assertEqual(self.client.post(url, {'type': 0, 'text': 'text annotation 3', 'start': 30}).status_code, 200)
        self.assertEqual(self.client.post(url, {'type': 1, 'text': 'text annotation 4', 'start': 40}).status_code, 200)
        self.assertEqual(self.client.post(url, {'type': 2, 'text': 'text annotation 5', 'start': 50}).status_code, 200)
        self.assertEqual(Annotation.objects.count(), 6)
        # Then submit the document for the Quest
        response = self.client.post(reverse('common:doc-quest-submit', kwargs={'quest_pk': self.task.pk, 'document_pk': doc.pk}), follow=True)
        self.client.logout()

        # As Anon user, export the documents submissions
        res = self.client.get(reverse('document:read-users-bioc', kwargs={'pubmed_id': doc.document_id, 'format_type': 'xml'}), follow=True)
        self.assertEqual(res.status_code, 200)
        bioc = BioCReader(source=res.content)
        bioc.read()

        # Make sure the BioC document has the opponent's infp
        self.assertEqual(len(bioc.collection.documents), 1)
        self.assertEqual(int(bioc.collection.documents[0].id), doc.document_id)
        self.assertEqual(len(bioc.collection.documents[0].passages), 2)
        self.assertEqual(len(bioc.collection.documents[0].passages[0].annotations), 0)
        self.assertEqual(len(bioc.collection.documents[0].passages[1].annotations), 6)
