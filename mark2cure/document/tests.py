from django.test import TestCase
from rest_framework.test import APIRequestFactory

from .tasks import get_pubmed_document
from .models import Document, Section, Pubtator

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

    def setUp(self):
        self.doc = Document.objects.first()

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
                document=self,
                session_id='',
                content__isnull=False)

        print 'QUERY: ', pub_query_set.count()

        response = self.client.get('/document/pubtator/{pmid}.json'.format(pmid=self.doc.document_id))
        json_string = response.content
        self.assertNotEqual(json_string, '', msg='API returned empty response for document BioC Pubtator Representation.')

        data = json.loads(json_string)

        # Make sure it's the same document in BioC as DB
        self.assertEqual(int(data.get('collection').get('document').get('id')), self.doc.document_id)
        self.assertEqual(len(data.get('collection').get('document').get('passage')), 2)
        self.assertEqual(data.get('collection').get('document').get('passage')[0].get('text'), self.doc.section_set.first().text)
        self.assertEqual(data.get('collection').get('document').get('passage')[1].get('text'), self.doc.section_set.last().text)

        # Make sure it doesn't contain any annotations
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

    def test_document_as_bioc_for_pairing(self):
        # (TODO) requires auth & previous submissions
        # Also requires GM testing
        # Could really just separate out all together
        pass

    def test_document_as_bioc_with_m2c(self):
        pass

