from django.core.urlresolvers import reverse
from django.test import TestCase

from ..test_base.test_base import TestBase
from ..document.models import Annotation
from ..common.models import Group
from ..task.models import Task

from .views import group_users_bioc

from ..common.bioc import BioCReader

import json


class GroupBioCAPIViews(TestCase, TestBase):
    # no templates (no json), directly to PA
    """class GroupBioCAPIViews can now inherit methods from TestBase,
    which is a class located in mark2cure/test_base/test_base.py
    """
    fixtures = ['tests_documents.json', 'tests_common.json']

    @classmethod
    def setUp(cls):

        cls.user_names = ['API-Test-User1', 'API-Test-User2', 'API-Test-User3',
        'API-Test-User4', 'API-Test-User5']

        cls.group = Group.objects.first()
        cls.task = Task.objects.first()

        # list of all the users we are creating
        cls.users = {}
        cls.user_annotation_list = []

    def test_quest_group_list(self):
        self.login_test_user('test_player')
        group = Group.objects.first()

        response = self.client.get(reverse('api:quest-group-api', kwargs={'group_pk': group.pk}))

        # Confirm API is online
        self.assertEqual(response.status_code, 200)

        # Confirm API provides Document PKs
        data = json.loads(response.content)
        data_string = str(data)
        # check that the exact pubmed IDs are in the json object
        self.assertTrue("1817, 1816, 1815" in data_string)
        self.client.logout()

    def test_group_users_bioc(self):
        self.create_new_user_accounts(self.user_names)
        # Get one document only (for users to share) (here just use user 0)
        # TODO: this could be improved. Currently using client to login and
        # retrieve a document to "share". Need to know how to do this more
        # elegently. JF 7/15/15
        doc = self.get_document_from_response(self.user_names[0])
        # Load fake annotations for all the users
        self.load_fake_annotations(self.user_names, doc)

        # Fetch the Group BioC as JSON to ensure is online
        response = self.client.get(reverse('api:group-users-bioc',
                                           kwargs={'group_pk': self.group.pk,
                                                   'format_type': 'json'}))
        self.assertEqual(response.status_code, 200)

        # Fetch the Group BioC for all user annotations
        response = self.client.get(reverse('api:group-users-bioc',
                                           kwargs={'group_pk': self.group.pk,
                                                   'format_type': 'xml'}))
        self.assertEqual(response.status_code, 200)
        r = BioCReader(source=response.content)
        r.read()

        # Does BioC have correct number of Group Documents
        self.assertEqual(len(r.collection.documents),
                         self.group.get_documents().count())

        # Does BioC have correct number of Group Annotations
        total_bioc_annotation_int = 0
        for bioc_doc in r.collection.documents:
            for bioc_passage in bioc_doc.passages:
                total_bioc_annotation_int += len(bioc_passage.annotations)
        self.assertEqual(Annotation.objects.count(), total_bioc_annotation_int)

    def test_group_pubtator_bioc(self):
        '''
        # As Anon user, export the documents submissions
        res = self.client.get(reverse('document:read-users-bioc',
                              kwargs={'pubmed_id': doc.document_id,
                              'format_type': 'xml'}), follow=True)
        self.assertEqual(res.status_code, 200)
        bioc = BioCReader(source=res.content)
        bioc.read()

        # Make sure the BioC document has the opponent's infp
        self.assertEqual(len(bioc.collection.documents), 1)
        self.assertEqual(int(bioc.collection.documents[0].id), doc.document_id)
        self.assertEqual(len(bioc.collection.documents[0].passages), 2)
        self.assertEqual(len(bioc.collection.documents[0].passages[0].annotations), 0)
        self.assertEqual(len(bioc.collection.documents[0].passages[1].annotations), 6)
        '''
        pass

    def test_group_list(self):
        # Confirm login required
        response = self.client.get(reverse('api:groups-api'))
        self.assertEqual(response.status_code, 302)

        self.login_test_user('test_player')
        response = self.client.get(reverse('api:groups-api'))

        # Confirm API is online
        self.assertEqual(response.status_code, 200)

        # Confirm API provides correct number of groups
        data = json.loads(response.content)
        self.assertEqual(len(data), Group.objects.count())

        self.client.logout()

    def test_leaderboard_users(self):
        self.login_test_user('test_player')

        for window in [1, 7, 30, 100]:
            # Confirm API is online
            response = self.client.get(reverse('api:leaderboard-users', kwargs={'day_window': window}))
            self.assertEqual(response.status_code, 200)

        self.client.logout()

    def test_leaderboard_teams(self):
        self.login_test_user('test_player')

        for window in [1, 7, 30, 100]:
            # Confirm API is online
            response = self.client.get(reverse('api:leaderboard-teams', kwargs={'day_window': window}))
            self.assertEqual(response.status_code, 200)

        self.client.logout()

class GroupUsersBioC(TestCase, TestBase):
    fixtures = ['large_group_3/data.json']

    def test_api_online(self):

        # Confirm API is online from direct access
        req = group_users_bioc({}, group_pk, 'xml')
        print req.content

        # Confirm API is online from URL
        response = self.client.get(reverse('api:quest-group-api',
            kwargs={'group_pk': 3}))

        self.assertEqual(response.status_code, 200)


