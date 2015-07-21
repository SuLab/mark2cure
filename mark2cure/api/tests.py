from django.test import TestCase
from django.core.urlresolvers import reverse

from ..document.models import Annotation
from ..common.models import Group, Task
from ..test_base.test_base import TestBase

from mark2cure.common.bioc import BioCReader


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

    def test_group_for_all_user_annotations(self):
        self.create_new_user_accounts(self.user_names, self.users)
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

    def test_group_for_all_pubtator_annotations(self):
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
