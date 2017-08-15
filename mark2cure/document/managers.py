from django.contrib.contenttypes.models import ContentType
from django.db import models, connection
from typing import List, Dict

from mark2cure.task.relation import relation_data_flat
from ..task.entity_recognition.models import EntityRecognitionAnnotation

import xml.etree.ElementTree as ET
from itertools import groupby
import pandas as pd

DF_COLUMNS = ('uid', 'source', 'user_id',
              'ann_type_idx', 'text',
              'document_pk', 'section_id', 'section_offset', 'offset_relative',
              'start_position', 'length')
APPROVED_TYPES = ['disease', 'gene_protein', 'drug']


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
                for pubtator_type_idx in range(len(pubtators[document_idx]) if document_idx < len(pubtators) else 0):
                    try:
                        root = ET.fromstring(pubtators[document_idx][pubtator_type_idx])
                        passage = root.findall(".//passage")[section_idx]
                        offset = int(passage.find("./offset").text)
                        for ann in passage.findall("./annotation"):
                            if 'Species' not in [infon.text for infon in ann.findall("./infon[@key='type']")]:
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

        return response

    def relation_df(self, documents=[], users=[]):
        ct = ContentType.objects.get(model='relationannotation')

        if len(documents):
            from .models import Document
            doc_arr = []
            for d in documents:
                if type(d) == Document:
                    doc_arr.append(str(d.pk))
                elif type(d) is str and d.isdigit():
                    doc_arr.append(d)
                elif type(d) is int:
                    doc_arr.append(str(d))
            filter_doc_level = 'WHERE `relation`.`document_id` IN ({0})'.format(','.join(doc_arr))
        else:
            filter_doc_level = ''

        if len(users):
            from django.contrib.auth.models import User
            user_arr = []
            for u in users:
                if type(u) == User:
                    user_arr.append(str(u.pk))
                elif type(u) is str and d.isdigit():
                    user_arr.append(u)
                elif type(u) is int:
                    user_arr.append(str(u))

            filter_user_level = '{0} `doc_view`.`user_id` IN ({1})'.format(
                'WHERE' if filter_doc_level == '' else 'AND',
                ','.join(user_arr))
        else:
            filter_user_level = ''

        cmd_str = ""
        with open('mark2cure/document/commands/get-relations-results.sql', 'r') as f:
            cmd_str = f.read()
        cmd_str = cmd_str.format(
            content_type_pk=ct.pk,
            filter_doc_level=filter_doc_level,
            filter_user_level=filter_user_level)

        c = connection.cursor()
        try:
            c.execute(cmd_str)

            REL_DF_COLUMNS = ('relation_id',
              'document_pk', 'document_pmid',
              'user_id',
              'relation_type', 'concept_1_id', 'concept_2_id',
              'answer', 'answer_text', 'created')

            df_arr = []
            for x in c.fetchall():
                df_arr.append({
                    'relation_id': x[0],
                    'document_pk': x[1], 'document_pmid': x[2],
                    'user_id': x[3],
                    'relation_type': x[4], 'concept_1_id': x[5], 'concept_2_id': x[6],
                    'answer': x[7], 'answer_text': filter(lambda d: d['id'] == x[7], relation_data_flat)[0]['text'], 'created': x[8]})

            return pd.DataFrame(df_arr, columns=REL_DF_COLUMNS)

        finally:
            c.close()

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

    def ner_df(self, documents=[], users=[], include_pubtator=True):
        if len(documents):
            from .models import Document
            doc_arr = []
            for d in documents:
                if type(d) == Document:
                    doc_arr.append(str(d.pk))
                elif type(d) is str and d.isdigit():
                    doc_arr.append(d)
                elif type(d) is int:
                    doc_arr.append(str(d))
            filter_doc_level = 'WHERE `document_section`.`document_id` IN ({0})'.format(','.join(doc_arr))
        else:
            filter_doc_level = ''

        if len(users):
            from django.contrib.auth.models import User
            user_arr = []
            for u in users:
                if type(u) == User:
                    user_arr.append(str(u.pk))
                elif type(u) is str and d.isdigit():
                    user_arr.append(u)
                elif type(u) is int:
                    user_arr.append(str(u))

            filter_user_level = '{0} `document_view`.`user_id` IN ({1})'.format(
                'WHERE' if filter_doc_level == '' else 'AND',
                ','.join(user_arr))
        else:
            filter_user_level = ''

        content_type_id = str(ContentType.objects.get_for_model(
            EntityRecognitionAnnotation.objects.first()).id)

        df_arr = []

        cmd_str = ""
        with open('mark2cure/document/commands/get-ner-results.sql', 'r') as f:
            cmd_str = f.read()
        cmd_str = cmd_str.format(content_type_pk=content_type_id, filter_doc_level=filter_doc_level, filter_user_level=filter_user_level)

        c = connection.cursor()
        try:
            c.execute(cmd_str)

            # Get the full writer in advnaced!!
            document_json = Document.objects.as_json(documents=documents)

            res = [x for x in c.fetchall()]

            # We group the response to reduce BioCDocument offset dict lookups
            for idx, doc_info in enumerate(groupby(res, lambda x: x[5])):
                document_pk, doc_group = doc_info

                # If a pubtator doesn't exist for the document, we can't include any annotations as the passage offsets need to come from Pubtator
                if document_pk == document_json[idx]['pk']:

                    # Use the BioC pubtator file for the offset values
                    offset_dict = {}
                    for passage in document_json[idx]['passages']:
                        offset_dict[int(passage['pk'])] = passage['offset']

                    # [(p['pk'], p['offset']) for p in res[0]['passages']]

                    for x in doc_group:
                        df_arr.append(self._create_er_df_row(
                            uid=x[0], source='db', user_id=x[8],
                            text=x[2], ann_type_idx=x[1],
                            document_pk=x[5], section_id=x[7], section_offset=offset_dict[x[7]], offset_relative=True,
                            start_position=x[3], length=len(x[2])))

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
            cmd_str = cmd_str.format(','.join(doc_arr))

            c = connection.cursor()
            try:
                c.execute(cmd_str)
                res = [x for x in c.fetchall()]
            finally:
                c.close()

            # Counter({'Disease': 3676, 'Chemical': 2928, 'Species': 1553, 'Gene': 1544, 'FamilyName': 536, 'DomainMotif': 20}) (Sampleing from DB 11/30/2016)
            pubtator_types = ['Disease', 'Gene', 'Chemical']
            for pubtator_content in res:
                r = BioCReader(source=pubtator_content[2])
                r.read()
                bioc_document = r.collection.documents[0]

                section_ids = pubtator_content[3].split(',')

                # Iterate over all the annotations in both passages
                for p_idx, passage in enumerate(bioc_document.passages):
                    for annotation in passage.annotations:

                        # Determine some meta-data (UID info) about the BioCAnnotation
                        annotation_type = None
                        uid_type = None
                        uid = None
                        for key in annotation.infons.keys():
                            if key == 'type':
                                annotation_type = annotation.infons.get(key, None)
                            else:
                                uid_type = key
                                uid = annotation.infons.get(uid_type, None)

                        # We're only interested in Pubtator Annotations that are the same concepts users highlight
                        if annotation_type in pubtator_types:
                            start, length = str(annotation.locations[0]).split(':')
                            df_arr.append(self._create_er_df_row(
                                uid=uid, source=uid_type if uid_type else None, user_id=None,
                                text=annotation.text, ann_type_idx=pubtator_types.index(annotation_type),
                                document_pk=pubtator_content[1], section_id=section_ids[p_idx], section_offset=passage.offset, offset_relative=False,
                                start_position=start, length=length))

        return pd.DataFrame(df_arr, columns=DF_COLUMNS)

