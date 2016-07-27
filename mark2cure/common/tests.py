from django.core.urlresolvers import reverse
from django.test import TestCase

from ..test_base.test_base import TestBase
from ..common.models import Group


class CommonViews(TestCase, TestBase):
    fixtures = ['tests_documents.json', 'tests_common.json']

    def test_home(self):
        # Confirm non-auth'd redirect
        response = self.client.get(reverse('common:home'))
        self.assertEqual(response.status_code, 200)
        # Ensure layout is the same
        self.assertInHTML('<meta name="application-name" content="Mark2Cure"/>', response.content)
        self.assertInHTML('<h2 class="modal-title">Subscribe</h2>', response.content)
        self.assertInHTML('<li><a href="http://sulab.org/category/mark2cure/">Blog</a>', response.content)

        # Confirm view online
        self.login_test_user('test_player')
        response = self.client.get(reverse('common:home'))
        self.assertEqual(response.status_code, 302)
        self.client.logout()

    def test_beta(self):
        # Confirm redirect
        response = self.client.get(reverse('common:beta'))
        self.assertEqual(response.status_code, 302)

    def test_why_i_mark2cure(self):
        # Confirm view online
        response = self.client.get(reverse('common:why-mark2cure'))
        self.assertEqual(response.status_code, 200)

        # Ensure footers appear on this webpage
        self.assert_footers_in_html(response.content)
        self.assertInHTML('<h2 class="text-xs-center">Why do you Mark2Cure?</h2>',
                          response.content)

    def test_dashboard(self):
        # Confirm non-auth'd redirect
        response = self.client.get(reverse('common:dashboard'))
        self.assertEqual(response.status_code, 302)

        # Confirm view online
        self.login_test_user('test_player')
        response = self.client.get(reverse('common:dashboard'))
        self.assertEqual(response.status_code, 200)
        # Check for standard footers and dashboard specific features.
        self.assert_footers_in_html(response.content)
        self.assertInHTML('<h4 id="myModalLabel" class="modal-title">Quest Instructions</h4>', response.content)
        self.assertInHTML('<h2 class="text-xs-center">Community Dashboard</h2>', response.content)
        self.client.logout()

    def test_group_view(self):
        group = Group.objects.first()

        # Confirm non-auth'd redirect
        response = self.client.get(reverse('common:group', kwargs={'group_stub': group.stub}))
        self.assertEqual(response.status_code, 302)

        # Confirm view online
        self.login_test_user('test_player')
        response = self.client.get(reverse('common:group', kwargs={'group_stub': group.stub}))
        self.assertEqual(response.status_code, 200)
        self.assert_footers_in_html(response.content)
        self.client.logout()

    def test_support(self):
        pass


class QuestViews(TestCase, TestBase):
    '''
        (TODO) This needs to be broken out into a
        'concept recognition' app and associated tests
    '''
    fixtures = ['tests_documents.json', 'tests_common.json']

    def test_quest_read_doc_results(self):
        pass

    def test_quest_read_doc_results_bioc(self):
        pass

    def test_quest_read_doc(self):
        pass

    def test_document_quest_submit(self):
        pass

    def test_quest_read(self):
        pass

    def test_quest_feedback(self):
        pass

    def test_quest_submit(self):
        pass
