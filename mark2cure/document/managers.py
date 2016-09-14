from django.db import models, connection
from django.contrib.contenttypes.models import ContentType
from ..task.entity_recognition.models import EntityRecognitionAnnotation
from mark2cure.task.relation import relation_data_flat
from ..common.bioc import BioCReader
from ..common.formatter import bioc_writer
import pandas as pd
from itertools import groupby


DF_COLUMNS = ('uid', 'source', 'user_id',
              'ann_type', 'text',
              'document_pk', 'section_id', 'section_offset', 'offset_relative',
              'start_position', 'length')
APPROVED_TYPES = ['disease', 'gene_protein', 'drug']


class DocumentManager(models.Manager):

    def as_writer(self, documents=[]):
        '''
            Return a blank BioC Writer that is based off the pubtator content.

            Problems: This requires every document to have at least 1 pubtator model
            Pros: This prevents us from generating our own BioC file which may
            have inconsistencies
        '''
        if len(documents):
            from .models import Document
            doc_arr = []
            for d in documents:
                if type(d) == Document:
                    doc_arr.append(str(d.pk))
                elif type(d) is str or type(d) is unicode and d.isdigit():
                    doc_arr.append(d)
                elif type(d) is int or type(d) is long:
                    doc_arr.append(str(d))

        cmd_str = '''
            SELECT
                `document_pubtator`.`document_id`,
                ANY_VALUE(`document_pubtator`.`content`),
                GROUP_CONCAT(DISTINCT `document_section`.`id`) as `section_ids`

            FROM `document_pubtator`

            JOIN `document_section`
                ON `document_section`.`document_id` = `document_pubtator`.`document_id`

            WHERE `document_pubtator`.`content` != ''
                AND `document_pubtator`.`session_id` = ''
                AND `document_pubtator`.`document_id` IN ({0})

            GROUP BY `document_pubtator`.`document_id`;
        '''.format(','.join(doc_arr))
        c = connection.cursor()
        try:
            c.execute(cmd_str)
            res = [(x[0], x[1], x[2]) for x in c.fetchall()]
        finally:
            c.close()

        writer = bioc_writer(None)
        for pubtator_content in res:
            section_ids = pubtator_content[2].split(',')
            r = BioCReader(source=pubtator_content[1])
            r.read()

            doc = r.collection.documents[0]
            doc.put_infon('document_pk', str(pubtator_content[0]))
            for idx, passage in enumerate(doc.passages):
                passage.clear_annotations()

                passage.put_infon('section', ['title', 'paragraph'][idx])
                passage.put_infon('id', str(section_ids[idx]))

            writer.collection.add_document(doc)

        return writer

    def relation_df(self, documents=[], users=[]):
        ct = ContentType.objects.get(model='relationannotation')

        if len(documents):
            from .models import Document
            doc_arr = []
            for d in documents:
                if type(d) == Document:
                    doc_arr.append(str(d.pk))
                elif type(d) is str or type(d) is unicode and d.isdigit():
                    doc_arr.append(d)
                elif type(d) is int or type(d) is long:
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
                elif type(u) is str or type(u) is unicode and d.isdigit():
                    user_arr.append(u)
                elif type(u) is int:
                    user_arr.append(str(u))

            filter_user_level = '{0} `doc_view`.`user_id` IN ({1})'.format(
                'WHERE' if filter_doc_level == '' else 'AND',
                ','.join(user_arr))
        else:
            filter_user_level = ''

        cmd_str = '''
            SELECT
                    `relation`.`id` as `relation_id`,
                    `doc_document`.`id` as `doc_pk`,
                    `doc_document`.`document_id` as `pmid`,
                    `doc_view`.`user_id`,

                    `relation`.`relation_type`,
                    `relation`.`concept_1_id`,
                    `relation`.`concept_2_id`,
                    `rel_anns`.`answer`,
                    `doc_ann`.`created`

            FROM `relation_relation` as `relation`

            INNER JOIN `relation_relationannotation` as `rel_anns`
                    ON `rel_anns`.`relation_id` = `relation`.`id`

            INNER JOIN `document_annotation` as `doc_ann`
                    ON `doc_ann`.`content_type_id` = {content_type_pk} AND `doc_ann`.`object_id` = `rel_anns`.`id`

            INNER JOIN `document_view` as `doc_view`
                    ON `doc_view`.`id` = `doc_ann`.`view_id`

            INNER JOIN `document_section` as `doc_section`
                    ON `doc_section`.`id` = `doc_view`.`section_id`

            INNER JOIN `document_document` as `doc_document`
                    ON `doc_document`.`id` = `doc_section`.`document_id`

            {filter_doc_level}
            {filter_user_level}

            ORDER BY `relation`.`id`
        '''.format(content_type_pk=ct.pk, filter_doc_level=filter_doc_level, filter_user_level=filter_user_level)
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

            '''
            for ro in df_arr:
                # (TODO) Select the longest text
                cdr_query = ConceptDocumentRelationship.objects.filter(document=relation.document)
                cdr1 = cdr_query.filter(concept_text__concept_id=relation.concept_1).first()
                cdr2 = cdr_query.filter(concept_text__concept_id=relation.concept_2).first()
            '''

            return pd.DataFrame(df_arr, columns=REL_DF_COLUMNS)

        finally:
            c.close()

    def _create_er_df_row(self,
                      uid, source='db', user_id=None,
                      ann_type='', text='',
                      document_pk=0, section_id=0, section_offset=0, offset_relative=True,
                      start_position=0, length=0):
        '''
            When offset_relative is False:
                start position is relative to the entire document and not the
                section it's contained within

            user_id can be None

            (TODO) Ann Types needs to be normalized
        '''
        ann_type = ann_type.lower()
        return {
            'uid': str(uid), 'source': str(source), 'user_id': int(user_id) if user_id else None,
            'ann_type': str(ann_type), 'text': str(text),
            'document_pk': int(document_pk), 'section_id': int(section_id), 'section_offset': int(section_offset), 'offset_relative': bool(offset_relative),
            'start_position': int(start_position), 'length': int(length)
        }

    def entity_recognition_df(self, documents=[], users=[], include_pubtator=True, writer=None):
        if len(documents):
            from .models import Document
            doc_arr = []
            for d in documents:
                if type(d) == Document:
                    doc_arr.append(str(d.pk))
                elif type(d) is str or type(d) is unicode and d.isdigit():
                    doc_arr.append(d)
                elif type(d) is int or type(d) is long:
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
                elif type(u) is str or type(u) is unicode and d.isdigit():
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

        cmd_str = '''
            SELECT
                `entity_recognition_entityrecognitionannotation`.`id`,
                `entity_recognition_entityrecognitionannotation`.`type`,
                `entity_recognition_entityrecognitionannotation`.`text`,
                `entity_recognition_entityrecognitionannotation`.`start`,
                `document_annotation`.`created`,
                `document_document`.`id` as `document_pk`,
                `document_document`.`document_id` as `pmid`,
                `document_view`.`section_id`,
                `document_view`.`user_id`

            FROM `entity_recognition_entityrecognitionannotation`

            INNER JOIN `document_annotation`
                ON `document_annotation`.`object_id` = `entity_recognition_entityrecognitionannotation`.`id` AND `document_annotation`.`content_type_id` = {content_type_pk}

            INNER JOIN `document_view`
                ON `document_annotation`.`view_id` = `document_view`.`id`

            INNER JOIN `document_section`
                ON `document_view`.`section_id` = `document_section`.`id`

            INNER JOIN `document_document`
                ON `document_document`.`id` = `document_section`.`document_id`

            {filter_doc_level}
            {filter_user_level}
        '''.format(content_type_pk=content_type_id, filter_doc_level=filter_doc_level, filter_user_level=filter_user_level)
        c = connection.cursor()
        try:
            c.execute(cmd_str)

            # Get the full writer in advnaced!!
            if not writer:
                writer = Document.objects.as_writer(documents=documents)

            res = [x for x in c.fetchall()]

            # We group the response to reduce BioCDocument offset dict lookups
            for key, doc_group in groupby(res, lambda x: x[5]):

                bioc_documents = filter(lambda d: d.infons.get('document_pk') == str(key), writer.collection.documents)
                # If a pubtator doesn't exist for the document, we can't include any annotations as the passage offsets need to come from Pubtator
                if len(bioc_documents) == 1:

                    # Use the BioC pubtator file for the offset values
                    offset_dict = {}
                    for passage in bioc_documents[0].passages:
                        offset_dict[int(passage.infons.get('id'))] = passage.offset

                    for x in doc_group:
                        df_arr.append(self._create_er_df_row(
                            uid=x[0], source='db', user_id=x[8],
                            text=x[2], ann_type=x[1],
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
            cmd_str = '''
                SELECT
                    `document_pubtator`.`id`,
                    `document_pubtator`.`document_id`,
                    `document_pubtator`.`content`,
                    GROUP_CONCAT(DISTINCT `document_section`.`id`) as `section_ids`

                FROM `document_pubtator`

                JOIN `document_section`
                    ON `document_section`.`document_id` = `document_pubtator`.`document_id`

                WHERE `document_pubtator`.`content` != ''
                    AND `document_pubtator`.`session_id` = ''
                    AND `document_pubtator`.`document_id` IN ({0})

                GROUP BY `document_pubtator`.`id`
            '''.format(','.join(doc_arr))
            c = connection.cursor()
            try:
                c.execute(cmd_str)
                res = [x for x in c.fetchall()]
            finally:
                c.close()

            for pubtator_content in res:
                r = BioCReader(source=pubtator_content[2])
                r.read()
                bioc_document = r.collection.documents[0]

                section_ids = pubtator_content[3].split(',')

                # Iterate over all the annotations in both passages
                for p_idx, passage in enumerate(bioc_document.passages):
                    for annotation in passage.annotations:

                        # Determine some meta-data (UID info) about the BioCAnnotation
                        uid = None
                        for key in annotation.infons.keys():
                            if key == 'type':
                                annotation_type = annotation.infons.get(key, None)
                            else:
                                uid_type = key
                                uid = annotation.infons.get(uid_type, None)

                        start, length = str(annotation.locations[0]).split(':')
                        df_arr.append(self._create_er_df_row(

                            uid=uid, source=uid_type if uid_type else None, user_id=None,
                            text=annotation.text, ann_type=annotation_type if annotation_type else None,
                            document_pk=pubtator_content[1], section_id=section_ids[p_idx], section_offset=passage.offset, offset_relative=False,
                            start_position=start, length=length))

        return pd.DataFrame(df_arr, columns=DF_COLUMNS)

