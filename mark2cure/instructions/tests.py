from django.core.urlresolvers import reverse
from django.test import TestCase

from ..test_base.test_base import TestBase


class InstructionViews(TestCase, TestBase):

    def test_disease_marking(self):
        # Confirm view online
        response = self.client.get(reverse('instructions:disease-marking'))
        self.assertEqual(response.status_code, 200)
        self.assert_footers_in_html(response.content)
        self.assertInHTML('<h2 class="text-center">Disease Concept Marking Rules</h2>',
                          response.content)
        self.assertInHTML('<input type="hidden" name="referral" value="/instructions/disease-marking/">',
                          response.content)

    def test_gene_marking(self):
        # Confirm view online
        response = self.client.get(reverse('instructions:gene-marking'))
        self.assertEqual(response.status_code, 200)
        self.assert_footers_in_html(response.content)
        self.assertInHTML('<h2 class="text-center">Gene Concept Marking Rules</h2>',
                          response.content)
        self.assertInHTML('<input type="hidden" name="referral" value="/instructions/gene-marking/">',
                          response.content)

    def test_treatment_marking(self):
        # Confirm view online
        response = self.client.get(reverse('instructions:treatment-marking'))
        self.assertEqual(response.status_code, 200)
        self.assert_footers_in_html(response.content)
        self.assertInHTML('<h2 class="text-center">Treatment Concept Marking Rules</h2>',
                          response.content)
        self.assertInHTML('<input type="hidden" name="referral" value="/instructions/treatment-marking/">',
                          response.content)
