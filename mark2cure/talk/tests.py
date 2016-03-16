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
        self.client.logout()

        # Confirm Moderator can views
        response = self.client.get(reverse('talk:home', kwargs={'pubmed_id': document.document_id}))
        self.assertEqual(response.status_code, 200)

        # Confirm view online for NON-MODERATOR
        # (TODO) Complete document

        self.client.logout()

    def test_annotation_search(self):
        # Confirm non-auth'd redirect
        response = self.client.get(reverse('talk:annotation-search'))
        self.assertEqual(response.status_code, 302)

        # Confirm view online
        self.login_test_user('new_test_player')
        response = self.client.get(reverse('talk:annotation-search'))
        self.assertEqual(response.status_code, 200)
        self.assert_footers_in_html(response.content)
        html_content_list = ['<th>Annotation Occurances</th>',
                             '<th>PMID</th>',
                             '<th>Document</th>',
                             '<h4 id="myModalLabel" class="modal-title">Quest Instructions</h4>']
        for item in html_content_list:
            self.assertInHTML(item, response.content)

    def test_recent_discussion(self):
        # Confirm non-auth'd redirect
        response = self.client.get(reverse('talk:recent'))
        self.assertEqual(response.status_code, 302)

        # Confirm view online
        self.login_test_user('new_test_player')
        response = self.client.get(reverse('talk:recent'))
        self.assertEqual(response.status_code, 200)
        self.assert_footers_in_html(response.content)
        html_content_list = ['<h2>Most Discussed Documents</h2>',
                             '<h2>Recent Annotations</h2>',
                             '<h2>Recent Comments</h2>',
                             '<p>Talk Pages</p>']
        for item in html_content_list:
            self.assertInHTML(item, response.content)
