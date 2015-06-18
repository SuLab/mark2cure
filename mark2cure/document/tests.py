from django.test import TestCase
from rest_framework.test import APIRequestFactory

from .tasks import get_pubmed_document
from .models import Document, Section, Pubtator


class DocumentProcessing(TestCase):
    def setUp(self):
        # random PMID from NCBI_corpus_testing
        self.pmid = 9467011
        get_pubmed_document(self.pmid)
        self.doc = Document.objects.get(document_id=self.pmid)

    def test_doc_init(self):
        self.assertEqual(Document.objects.count(), 1)
        self.assertEqual(self.doc.document_id, self.pmid)

    def test_section_creation(self):
        sections = Section.objects.filter(document=self.doc)
        self.assertEqual(sections.count(), 2)

    def test_pubtator_init(self):
        pubtators = Pubtator.objects.filter(document=self.doc)
        self.assertEqual(pubtators.count(), 3)

    def test_doc_methods(self):
        self.assertEqual(self.doc.available_sections().count(), 2)
        self.assertEqual(self.doc.count_available_sections(), 2)

    def test_padding_processing(self):
        # Actually determine if a doc would be repadded
        self.assertEqual(self.doc.update_padding(), False)


