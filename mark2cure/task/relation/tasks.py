from .models import Concept, ConceptText, ConceptDocumentRelationship, Relation, RelationGroup
from ...document.models import Document
from ...common.formatter import clean_df

import itertools
from celery import task

"""
# makes concept lists for only documents that have at least two dictionaries
that contain concepts. Must be done server side so that it can be
determined if there are non-overlapping concepts
(at least one pair, or "relation pair" per document).
"""

stype_arr = ['d', 'g', 'c']


@task()
def import_concepts():
    """
        This is where the number of documents to prepopulate the relation app starts
        from. For each document we make concept dictionaries from each pubtator
            1) check for overlapping words
            2) filter out bad documents
            3) Make sure there are at least two concept_types so that pairs can be formed
            4) Do not use concepts without a unique identification (UID)
            5) Use the longest word as the representative "text" if there are multiple texts
            per UID
    """

    for group in RelationGroup.objects.all():

        for document in group.documents.all():

            df = Document.objects.entity_recognition_df(documents=[document], include_pubtator=True)
            # (TODO) Do we need protection or not? I'm confused.
            df = clean_df(df, overlap_protection=False, allow_duplicates=False)

            # We're only interested in Pubtator identifiers
            df = df[df['user_id'] == -1]

            # Put these into the DB
            for index, row in df.iterrows():
                # (TODO) Maybe make these bulk creates
                c, created = Concept.objects.get_or_create(id=row['uid'])
                ct, created = ConceptText.objects.get_or_create(concept_id=row['uid'], text=row['text'])
                # (TODO) Consider ann_type in addition to source
                # agree with this comment. 'source' is too variable. Just use ann_type -JF
                cdr, created = ConceptDocumentRelationship.objects.get_or_create(concept_text=ct, document=document, stype=stype_arr[row['ann_type_idx']])


@task()
def compute_relationships():
    """
        This method takes a document and a relation pair list and makes the
        appropriate "relation." The reason relations are stored instead of
        created on the fly, is because there was a significant amount of code
        performed to determine whether or not relations even exist.  We need
        those relation pairs to determine if the document is worth showing to
        a user.
    """

    approved_relationship_types = ['c_d', 'g_d', 'g_c']

    document_ids = set(ConceptDocumentRelationship.objects.values_list('document_id', flat=True))
    for document_id in document_ids:
        cdr_query = ConceptDocumentRelationship.objects.filter(document_id=document_id, stype__isnull=False)

        document_concepts_set = set(cdr_query.values_list('concept_text__concept_id', flat=True))

        uid_stype_hash = cdr_query.values('concept_text__concept_id', 'concept_text', 'stype')

        for uid_a, uid_b in itertools.combinations(document_concepts_set, 2):
            uid_a_dict = filter(lambda cdr: cdr['concept_text__concept_id'] == uid_a, uid_stype_hash)[0]
            uid_b_dict = filter(lambda cdr: cdr['concept_text__concept_id'] == uid_b, uid_stype_hash)[0]
            relation_type = uid_a_dict['stype'] + '_' + uid_b_dict['stype']

            # Keep Gene to the LEFT of Diseases
            # Keep Disease to the RIGHT (this is important for future analysis)
            if relation_type in [x[::-1] for x in approved_relationship_types]:
                uid_a, uid_b = uid_b, uid_a
                relation_type = relation_type[::-1]

            if relation_type in approved_relationship_types:
                Relation.objects.get_or_create(
                    document_id=document_id,
                    relation_type=relation_type,
                    concept_1_id=uid_a,
                    concept_2_id=uid_b)
