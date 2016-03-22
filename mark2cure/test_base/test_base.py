from django.core.urlresolvers import reverse
from django.contrib.auth.models import User

from ..document.models import Annotation

from brabeion import badges
from random import randint
import random


class TestBase(object):
    """ TestBase has methods that can be inherited from for unit testing"""

    def __init__(self):
        pass

    def load_fake_annotations(self, user_names, doc):
        """Create simulated user annotations
        Input is a document object from response.context, & list of user names
        """
        self.assertEqual(Annotation.objects.count(), 0)
        total_ann_count = 0

        # update id_generator to contain realistic biomedical text annotations.
        # JF 7/13/15
        def _id_generator(doc_text):
            split_doc_text = doc_text.split(" ")
            word_length = len(split_doc_text)
            user_annotation = []
            while user_annotation == []:
                # 3 selected for max of annotation text
                random_number_dif = randint(0, 3)
                random_number_end = randint(0, word_length)
                random_number_start = random_number_end - random_number_dif
                user_annotation = split_doc_text[random_number_start:random_number_end]
            user_annotation = " ".join(user_annotation)
            # return one annotation, can be added to list of annotations
            return user_annotation

        for user_name in user_names:
            # "Print" comments below for useful text annotation info
            # print user_name
            # print doc
            # Each user gets a document to annotate
            user_annotation_list = []
            for section in doc.available_sections():
                ann_count = random.randint(0, 30)
                for x in range(ann_count):
                    url = reverse('task-entity-recognition:create',
                                  kwargs={'task_pk': self.task.pk,
                                  'section_pk': section.pk})
                    # user_annotations shall be from the current section text
                    user_annotation = _id_generator(section.text)
                    self.assertEqual(self.client.post(url,
                                     {'type': random.randint(0, 2),
                                     'text': user_annotation,
                                     'start': random.randint(0, len(section.text))}).status_code, 200)
                    user_annotation_list.append(user_annotation)
                # print self.user_annotation_list
                total_ann_count = total_ann_count + ann_count
                self.assertEqual(Annotation.objects.count(), total_ann_count)

    def create_new_user_accounts(self, user_names):
        """ Input is a list of users to initiate first login"""
        for user_name in user_names:
            user_name = User.objects.create_user(user_name, password='password')
            badges.possibly_award_badge("skill_awarded",
                                        user=user_name,
                                        level=7, force=True)

    def login_test_user(self, users_name):
        self.create_new_user_accounts([users_name, ])
        user = User.objects.first()
        password = 'password'
        self.client.login(username=user.username, password=password)

    def get_document_from_response(self, user_name):
        """ If you need one document, arbitrarily use doc from user_name"""
        self.client.login(username=user_name, password='password')
        url = reverse('task-entity-recognition:quest-home', kwargs={'quest_pk': self.task.pk})
        response = self.client.get(url, follow=True)
        doc = response.context['document']
        return doc

    def assert_footers_in_html(self, response_content):
        footers = ['<p>About Mark2Cure</p>',
                   '<p>FAQ</p>',
                   '<p>Help / Contact</p>',
                   '<p>Terms of Service / Privacy Policy</p>']
        for footer in footers:
            self.assertInHTML(footer, response_content)
