from django.core.urlresolvers import reverse
from django.test import TestCase

from ..test_base.test_base import TestBase
from ..common.models import Group


class CommonViews(TestCase, TestBase):
    fixtures = ['tests_documents.json', 'tests_common.json']

    @classmethod
    def setUp(cls):
        # (JK TODO) remove this required global for TestBase inheritance
        cls.users = {}

    def test_home(self):
        # Confirm non-auth'd redirect
        response = self.client.get(reverse('common:home'))
        self.assertEqual(response.status_code, 200)

        # Confirm view online
        user, password = self.get_test_user()
        self.client.login(username=user.username, password=password)
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

    def test_dashboard(self):
        # Confirm non-auth'd redirect
        response = self.client.get(reverse('common:dashboard'))
        self.assertEqual(response.status_code, 302)

        # Confirm view online
        user, password = self.get_test_user()
        self.client.login(username=user.username, password=password)
        response = self.client.get(reverse('common:dashboard'))
        self.assertEqual(response.status_code, 200)
        self.client.logout()

    def test_group_view(self):
        group = Group.objects.first()

        # Confirm non-auth'd redirect
        response = self.client.get(reverse('common:group', kwargs={'group_stub': group.stub}))
        self.assertEqual(response.status_code, 302)

        # Confirm view online
        user, password = self.get_test_user()
        self.client.login(username=user.username, password=password)
        response = self.client.get(reverse('common:group', kwargs={'group_stub': group.stub}))
        self.assertEqual(response.status_code, 200)
        self.client.logout()

    def test_support(self):
        pass


class QuestViews(TestCase, TestBase):
    '''
        (TODO) This needs to be broken out into a
        'concept recognition' app and associated tests
    '''
    fixtures = ['tests_documents.json', 'tests_common.json']

    @classmethod
    def setUp(cls):
        cls.users = {}

    def test_quest_read_doc_results(self):
        pass

    def test_quest_read_doc_results_bioc(self):
        pass

    def quest_read_doc(self):
        pass

    def test_document_quest_submit(self):
        pass

    def test_quest_read(self):
        pass

    def test_quest_feedback(self):
        pass

    def test_quest_submit(self):
        pass

