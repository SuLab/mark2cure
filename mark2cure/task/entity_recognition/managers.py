from django.db import models


class EntityRecognitionAnnotationManager(models.Manager):

    def document_pks_by_text_and_document_pks(self, text, document_pks):
        res = self.raw("""
            SELECT  entity_recognition_entityrecognitionannotation.id,
                    document_section.document_id as document_id
            FROM entity_recognition_entityrecognitionannotation
            LEFT OUTER JOIN document_annotation
                ON document_annotation.object_id = entity_recognition_entityrecognitionannotation.id
                LEFT OUTER JOIN document_view
                    ON document_view.id = document_annotation.view_id
                        LEFT OUTER JOIN document_section
                            ON document_section.id = document_view.section_id
            WHERE ( document_section.document_id IN ({0})
                    AND entity_recognition_entityrecognitionannotation.text != ''
                    AND entity_recognition_entityrecognitionannotation.text = '{1}')
        """.format( ', '.join('\'' + str(item) + '\'' for item in document_pks) , text));
        return [x.document_id for x in res]


    def annotations_texts_by_created_and_document_pks(self, created_datetime, document_pks):
        res = self.raw("""
            SELECT  entity_recognition_entityrecognitionannotation.id,
                    entity_recognition_entityrecognitionannotation.text
            FROM entity_recognition_entityrecognitionannotation
            LEFT OUTER JOIN document_annotation
                ON document_annotation.object_id = entity_recognition_entityrecognitionannotation.id
                LEFT OUTER JOIN document_view
                    ON document_view.id = document_annotation.view_id
                        LEFT OUTER JOIN document_section
                            ON document_section.id = document_view.section_id
            WHERE ( document_section.document_id IN ({0})
                    AND entity_recognition_entityrecognitionannotation.text != ''
                    AND document_annotation.created >= '{1}')
        """.format( ', '.join('\'' + str(item) + '\'' for item in document_pks) , created_datetime));
        return [x.text for x in res]


    def annotations_texts_by_created(self, created_datetime):
        res = self.raw("""
            SELECT  entity_recognition_entityrecognitionannotation.id,
                    entity_recognition_entityrecognitionannotation.text
            FROM entity_recognition_entityrecognitionannotation
            LEFT OUTER JOIN document_annotation
                ON document_annotation.object_id = entity_recognition_entityrecognitionannotation.id
            WHERE ( entity_recognition_entityrecognitionannotation.text != ''
                    AND document_annotation.created >= '{0}')
        """.format(created_datetime));
        return [x.text for x in res]


    def annotations_texts_for_document_and_type(self, doc_pk, type):
        res = self.raw("""
            SELECT  entity_recognition_entityrecognitionannotation.id,
                    entity_recognition_entityrecognitionannotation.text
            FROM entity_recognition_entityrecognitionannotation
            LEFT OUTER JOIN document_annotation
                ON document_annotation.object_id = entity_recognition_entityrecognitionannotation.id
                LEFT OUTER JOIN document_view
                    ON document_view.id = document_annotation.view_id
                        LEFT OUTER JOIN document_section
                            ON document_section.id = document_view.section_id
                            LEFT OUTER JOIN document_document
                                ON document_document.id = document_section.document_id
            WHERE (document_document.id = {0} AND entity_recognition_entityrecognitionannotation.type = '{1}' AND entity_recognition_entityrecognitionannotation.text != '')
        """.format(doc_pk, type));
        return [x.text for x in res]

    def annotations_for_section_pk(self, section_pk):
        res = self.raw("""
            SELECT
                entity_recognition_entityrecognitionannotation.id,
                entity_recognition_entityrecognitionannotation.type,
                entity_recognition_entityrecognitionannotation.text,
                entity_recognition_entityrecognitionannotation.start,
                document_annotation.created,
                document_view.section_id,
                document_view.user_id
            FROM entity_recognition_entityrecognitionannotation
            LEFT OUTER JOIN document_annotation
                ON document_annotation.object_id = entity_recognition_entityrecognitionannotation.id
                LEFT OUTER JOIN document_view
                    ON document_annotation.view_id = document_view.id
            WHERE section_id = {0}
        """.format(section_pk));
        return res

    def annotations_for_view_pks(self, view_pks):
        res = self.raw("""
            SELECT  entity_recognition_entityrecognitionannotation.id,
                    entity_recognition_entityrecognitionannotation.text,
                    entity_recognition_entityrecognitionannotation.type,
                    entity_recognition_entityrecognitionannotation.start
            FROM entity_recognition_entityrecognitionannotation
            LEFT OUTER JOIN document_annotation
                ON document_annotation.object_id = entity_recognition_entityrecognitionannotation.id
                LEFT OUTER JOIN document_view
                    ON document_view.id = document_annotation.view_id
            WHERE document_view.id IN ({0})
        """.format( ', '.join('\'' + str(item) + '\'' for item in view_pks) ));
        return [x for x in res]

    def annotations_for_document_pmid(self, doc_pmid):
        res = self.raw("""
            SELECT
                entity_recognition_entityrecognitionannotation.id,
                entity_recognition_entityrecognitionannotation.type,
                entity_recognition_entityrecognitionannotation.text,
                entity_recognition_entityrecognitionannotation.start,
                document_annotation.created,
                document_view.section_id,
                document_view.user_id,
                document_document.document_id
            FROM entity_recognition_entityrecognitionannotation
            LEFT OUTER JOIN document_annotation
                ON document_annotation.object_id = entity_recognition_entityrecognitionannotation.id
                LEFT OUTER JOIN document_view
                    ON document_view.id = document_annotation.view_id
                        LEFT OUTER JOIN document_section
                            ON document_section.id = document_view.section_id
                            LEFT OUTER JOIN document_document
                                ON document_document.id = document_section.document_id
            WHERE document_document.document_id = {0}
        """.format( doc_pmid ));

        return [i for i in res]


    def annotations_for_section_and_user(self, section_pk, user_pk):
        res = self.raw("""
            SELECT
                entity_recognition_entityrecognitionannotation.id,
                entity_recognition_entityrecognitionannotation.type,
                entity_recognition_entityrecognitionannotation.text,
                entity_recognition_entityrecognitionannotation.start,
                document_annotation.created,
                document_view.section_id,
                document_view.user_id
            FROM entity_recognition_entityrecognitionannotation
            LEFT OUTER JOIN document_annotation
                ON document_annotation.object_id = entity_recognition_entityrecognitionannotation.id
                LEFT OUTER JOIN document_view
                    ON document_annotation.view_id = document_view.id
            WHERE (section_id = {0} AND user_id = {1})
        """.format(section_pk, user_pk));
        return res

