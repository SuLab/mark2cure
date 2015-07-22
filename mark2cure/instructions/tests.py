from django.core.urlresolvers import reverse
from django.test import TestCase

from ..test_base.test_base import TestBase


class InstructionViews(TestCase, TestBase):

    def test_disease_marking(self):
        # Confirm view online
        response = self.client.get( reverse('instructions:disease-marking') )
        self.assertEqual(response.status_code, 200)

    def test_gene_marking(self):
        # Confirm view online
        response = self.client.get( reverse('instructions:gene-marking') )
        self.assertEqual(response.status_code, 200)

    def test_treatment_marking(self):
        # Confirm view online
        response = self.client.get( reverse('instructions:treatment-marking') )
        self.assertEqual(response.status_code, 200)

