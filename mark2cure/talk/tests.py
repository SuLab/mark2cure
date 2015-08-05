from django.core.urlresolvers import reverse
from django.test import TestCase

from ..test_base.test_base import TestBase
from ..document.models import Document


class TalkViews(TestCase, TestBase):

    fixtures = ['tests_documents.json', 'tests_common.json']

    def test_home(self):
        document = Document.objects.first()

        # Confirm non-auth'd redirect
        response = self.client.get(reverse('talk:home', kwargs={'pubmed_id': document.document_id}))
        self.assertEqual(response.status_code, 302)

        # Confirm view decorator protects non complete users
        self.login_test_user('new_test_player')
        response = self.client.get(reverse('talk:home', kwargs={'pubmed_id': document.document_id}))
        self.assertEqual(response.status_code, 403)

        # Confirm view online
        # (TODO)

        self.client.logout()

    def test_annotation_search(self):
        # Confirm non-auth'd redirect
        response = self.client.get(reverse('talk:annotation-search'))
        self.assertEqual(response.status_code, 302)

        # Confirm view online
        self.login_test_user('new_test_player')
        response = self.client.get(reverse('talk:annotation-search'))
        self.assertEqual(response.status_code, 200)

    def test_recent_discussion(self):
        # Confirm non-auth'd redirect
        response = self.client.get(reverse('talk:recent'))
        self.assertEqual(response.status_code, 302)

        # Confirm view online
        self.login_test_user('new_test_player')
        response = self.client.get(reverse('talk:recent'))
        self.assertEqual(response.status_code, 200)
