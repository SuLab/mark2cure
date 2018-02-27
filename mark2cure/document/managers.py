from django.contrib.contenttypes.models import ContentType
from django.db import models, connection
from typing import List, Dict

from ..task.ner.models import EntityRecognitionAnnotation

import xml.etree.ElementTree as ET
from itertools import groupby
import pandas as pd

NER_DF_COLUMNS = ('uid', 'source', 'user_id',
                  'ann_type_idx', 'text',
                  'document_pk', 'section_id', 'section_offset', 'offset_relative',
                  'start_position', 'length')

RE_DF_COLUMNS = ('relation_id',
  'document_pk', 'document_pmid',
  'user_id',
  'relation_type', 'concept_1_id', 'concept_2_id',
  'answer', 'answer_text', 'created')

APPROVED_TYPES = ['disease', 'gene_protein', 'drug']

PUBTATOR_TYPES = ['Disease', 'Gene', 'Chemical']


class DocumentManager(models.Manager):

    def as_json(self, document_pks: List[int], pubtators=[]) -> List[Dict]:
        """Represent the selection of documents as a Array of (dict)Documents

        Args:
            documents_pks (list): The selection of JSON Documents to return
            pubtators (list): The paired pubtator content bodies

        Returns:
            list: The list of (dict)Documents
        """
        assert len(document_pks) >= 1, "No documents supplied to generator JSON"
        # assert len(pubtators) >= 1 and len(document_pks) != len(pubtators), "Incorrect pairing of Documents and pubtators."

        cmd_str = ""
        with open('mark2cure/document/commands/get-documents.sql', 'r') as f:
            cmd_str = f.read()
        cmd_str = cmd_str.format(','.join([str(x) for x in document_pks]))

        c = connection.cursor()
        try:
            c.execute(cmd_str)
            doc_queryset = [dict(zip(['pk', 'pmid', 'section',
                              'section_pk', 'text'], x)) for x in c.fetchall()]
        finally:
            c.close()

        response = []
        for document_idx, document_group in enumerate(groupby(doc_queryset, lambda x: x['pk'])):
            document_pk, document_sections = document_group

            passages = []
            for section_idx, section_dict in enumerate(document_sections):
                offset = 0

                passage_annotations = []
                # If pubtators are available, iterate over them
                for pubtator_type_idx in range(len(pubtators[document_idx]) if document_idx < len(pubtators) else 0):
                    try:
                        root = ET.fromstring(pubtators[document_idx][pubtator_type_idx])
                        passage = root.findall(".//passage")[section_idx]
                        offset = int(passage.find("./offset").text)
                        for ann in passage.findall("./annotation"):
                            if [infon.text for infon in ann.findall("./infon[@key='type']")][0] in PUBTATOR_TYPES:
                                passage_annotations.append({
                                    'type_id': pubtator_type_idx,
                                    'start': int(ann.find('location').attrib['offset']),
                                    'text': ann.find('text').text
                                })
                    except SyntaxError:
                        pass

                passages.append({
                    'section': section_dict['section'],
                    'pk': section_dict['section_pk'],
                    'text': section_dict['text'],
                    'offset': offset,
                    'annotations': passage_annotations
                })

            response.append({
                'pk': document_pk,
                'passages': passages
            })

        # If no pubtators, gotta fill in the offset ourselves
        for document_idx, document in enumerate(response):
            if len(pubtators) == 0:
                for passage_idx, passage in enumerate(document['passages']):
                    if passage_idx == 0:
                        continue
                    l_passage = document['passages'][passage_idx - 1]
                    passage['offset'] = l_passage['offset'] + len(l_passage['text'])

        return response

    def re_df(self, document_pks: List[int], user_pks: List[int]=[]):
        """Relationship Extraction Results DataFrame

        (TODO) Not finished

        Args:
            documents_pks (list): The selection of JSON Documents to return
            user_pks (list): The selection of Users to include in the results

        Returns:
            pd.DataFrame: The list of (dict)Documents
        """
        assert len(document_pks) >= 1, "No documents supplied to Relationship Extraction Dataframe"
        filter_doc_level = 'WHERE `relation`.`document_id` IN ({0})'.format(','.join([str(x) for x in document_pks]))

        if user_pks and len(user_pks):
            filter_user_level = '{0} `doc_view`.`user_id` IN ({1})'.format(
                'WHERE' if filter_doc_level == '' else 'AND',
                ','.join([str(x) for x in user_pks]))
        else:
            filter_user_level = ""

        cmd_str = ""
        with open('mark2cure/document/commands/get-relations-results.sql', 'r') as f:
            cmd_str = f.read()
        cmd_str = cmd_str.format(
            content_type_pk=ContentType.objects.get(model='relationannotation').pk,
            filter_doc_level=filter_doc_level,
            filter_user_level=filter_user_level)

        c = connection.cursor()
        try:
            c.execute(cmd_str)
            re_queryset = [dict(zip(['relation_id', 'document_pk', 'document_pmid',
                              'user_id', 'relation_type', 'concept_1_id',
                              'concept_2_id', 'answer', 'created'], x)) for x in c.fetchall()]

            print(re_queryset[0])

        # 'answer_text': filter(lambda d: d['id'] == x[7], relation_data_flat)[0]['text'],

        finally:
            c.close()

        df_arr = []
        return pd.DataFrame(df_arr, columns=RE_DF_COLUMNS)

    def _create_er_df_row(self,
                      uid, source='db', user_id=None,
                      ann_type_idx=0, text='',
                      document_pk=0, section_id=0, section_offset=0, offset_relative=True,
                      start_position=0, length=0):
        '''
            When offset_relative is False:
                start position is relative to the entire document and not the
                section it's contained within

            user_id can be None

        '''
        return {
            'uid': str(uid), 'source': str(source), 'user_id': int(user_id) if user_id else None,
            'ann_type_idx': int(ann_type_idx), 'text': str(text),
            'document_pk': int(document_pk), 'section_id': int(section_id), 'section_offset': int(section_offset), 'offset_relative': bool(offset_relative),
            'start_position': int(start_position), 'length': int(length)
        }

    def ner_df(self, document_pks: List[int], user_pks: List[int]=[], include_pubtator=True):
        """
        Returns:
            pd.DataFrame: Named Entity Recognition Results Dataframe
        """

        assert len(document_pks) >= 1, "No documents supplied to Relationship Extraction Dataframe"
        filter_doc_level = 'WHERE `document_section`.`document_id` IN ({0})'.format(','.join([str(x) for x in document_pks]))

        if len(user_pks):
            filter_user_level = '{0} `document_view`.`user_id` IN ({1})'.format(
                'WHERE' if filter_doc_level == '' else 'AND',
                ','.join([str(x) for x in user_pks]))
        else:
            filter_user_level = ''

        cmd_str = ""
        with open('mark2cure/document/commands/get-ner-results.sql', 'r') as f:
            cmd_str = f.read()
        cmd_str = cmd_str.format(
            content_type_pk=ContentType.objects.get_for_model(EntityRecognitionAnnotation).id,
            filter_doc_level=filter_doc_level,
            filter_user_level=filter_user_level)

        df_arr = []
        c = connection.cursor()
        try:
            c.execute(cmd_str)

            # Get the full writer in advnaced!!
            document_json = self.as_json(document_pks=document_pks)

            ner_queryset = [dict(zip(['pk', 'type_idx', 'text',
                              'start', 'created', 'document_pk',
                              'pmid', 'section_pk', 'user_id'], x)) for x in c.fetchall()]

            # We group the response to reduce offset dict lookups
            for document_idx, document_group in enumerate(groupby(ner_queryset, lambda x: x['document_pk'])):
                document_pk, document_annotations = document_group

                # If a pubtator doesn't exist for the document, we can't include any annotations as the passage offsets need to come from Pubtator
                if document_pk == document_json[document_idx]['pk']:

                    # Use the (dict)Document JSON file for the offset values
                    offset_dict = {}
                    for passage in document_json[document_idx]['passages']:
                        offset_dict[int(passage['pk'])] = passage['offset']

                    for annotation in document_annotations:
                        df_arr.append(self._create_er_df_row(
                            uid=annotation['pk'], source='db', user_id=annotation['user_id'],
                            text=annotation['text'], ann_type_idx=annotation['type_idx'],
                            document_pk=annotation['document_pk'], section_id=annotation['section_pk'], section_offset=offset_dict[annotation['section_pk']], offset_relative=True,
                            start_position=annotation['start'], length=len(annotation['text'])))

        finally:
            c.close()

        if include_pubtator:
            '''
                This is the component that merges the 3 different pubtator
                reponses into 1 main file. It performances selective
                ordering and precedence for some annotations types / instances
            '''
            cmd_str = ""
            with open('mark2cure/document/commands/get-ner-pubtator-results.sql', 'r') as f:
                cmd_str = f.read()
            cmd_str = cmd_str.format(','.join([str(x) for x in document_pks]))

            c = connection.cursor()
            try:
                c.execute(cmd_str)
                pubtator_queryset = [dict(zip(['pk', 'document_pk', 'content',
                                  'section_pks'], x)) for x in c.fetchall()]
            finally:
                c.close()

            for pubtator in pubtator_queryset:
                section_ids = pubtator['section_pks'].split(',')
                root = ET.fromstring(pubtator['content'])

                # bioc_document = r.collection.documents[0]

                # Iterate over all the annotations in both passages
                for passage_idx, passage in enumerate(root.findall(".//passage")):
                    offset = int(passage.find("./offset").text)

                    for annotation in passage.findall("./annotation"):
                        annotation_type = [infon.text for infon in annotation.findall("./infon[@key='type']")][0]
                        uids_list = [infon.text for infon in annotation.findall("./infon[@key='identifier']")]
                        uid = uids_list[0] if uids_list else None

                        # We're only interested in Pubtator Annotations that are the same concepts users highlight
                        if annotation_type in PUBTATOR_TYPES:
                            start = int(annotation.find('location').attrib['offset'])
                            text = annotation.find('text').text

                            df_arr.append(self._create_er_df_row(
                                uid=uid, source='identifier' if uid else None, user_id=None,
                                text=text, ann_type_idx=PUBTATOR_TYPES.index(annotation_type),
                                document_pk=pubtator['document_pk'], section_id=section_ids[passage_idx], section_offset=offset, offset_relative=False,
                                start_position=start, length=len(text)))

        return pd.DataFrame(df_arr, columns=NER_DF_COLUMNS)

