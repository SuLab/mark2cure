from django.test import TestCase
from rest_framework.test import APIRequestFactory
from django.core.urlresolvers import reverse

from ..document.tasks import get_pubmed_document
from ..document.models import Document, Section, Pubtator, Annotation, View
from ..common.models import Group, Task, UserQuestRelationship
from django.contrib.auth.models import User
from brabeion import badges
##from ..api.tests import GroupBioCAPIViews ### Jennifer TODO

from mark2cure.common.bioc import BioCReader
from Bio import Entrez, Medline
import datetime
import random
from random import randint
import json


class TestBase(object):
    """ TestBase has methods that can be inherited from for unit testing
    """

    def __init__(self):
        pass

    def load_fake_annotations(self):

        self.assertEqual(Annotation.objects.count(), 0)
        total_ann_count = 0

        # update id_generator to contain realistic biomedical text annotations. JF 7/13/15
        def _id_generator(doc_text):
            split_doc_text = doc_text.split(" ")
            word_length = len(split_doc_text)
            user_annotation = []
            while user_annotation == []:
                random_number_dif = randint(0,3) # 3 selected for max of annotation text
                random_number_end = randint(0,word_length)
                random_number_start = random_number_end - random_number_dif
                user_annotation = split_doc_text[random_number_start:random_number_end]
            user_annotation = " ".join(user_annotation)
            return user_annotation

        for user_name in self.user_names:
            # Jennifer TODO
            # Submit Annotations (As User 1) so they show up when inspecting the M2C submissions
            self.client.login(username=user_name, password='password')
            response = self.client.get(reverse('common:quest-home', kwargs={'quest_pk': self.task.pk}), follow=True)

            # Each User will get a different document from the group to annotate
            doc = response.context['document']
            self.user_annotation_list = []
            for section in doc.available_sections():
                # number of user annotations can range from 0 to 30
                ann_count = random.randint(0,30)
                for x in range(ann_count):
                    url = reverse('document:create', kwargs={'task_pk': self.task.pk, 'section_pk': section.pk})
                    # user_annotations from the current section text
                    user_annotation = _id_generator(section.text)
                    # TODO: Max, is this itself an assert statemn to check if the url is working?
                    self.assertEqual(self.client.post(url, {'type': random.randint(0,2), 'text': user_annotation, 'start': random.randint(0, len(section.text))}).status_code, 200)
                    self.user_annotation_list.append(user_annotation)
                total_ann_count = total_ann_count + ann_count
                self.assertEqual(Annotation.objects.count(), total_ann_count)
            # prints final annotation list for each user for testing
            print self.user_annotation_list
