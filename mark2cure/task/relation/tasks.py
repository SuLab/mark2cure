from ...common.bioc import BioCReader
from ...document.models import Document, Pubtator
from .models import Concept, ConceptText, ConceptDocumentRelationship, Relation

import re
import itertools
from celery import task
import pandas as pd

"""
# makes concept lists for only documents that have at least two dictionaries
that contain concepts. Must be done server side so that it can be
determined if there are non-overlapping concepts
(at least one pair, or "relation pair" per document).
"""

stype_hash = {
    'MESH': 'c',
    'MEDIC': 'd',
    'NCBI Gene': 'g'
}

def check_for_overlaps(c1_dict_locations, c2_dict_locations):
    """ Takes the location from the concepts and checks to see
    if they overlap. Swaps c1 and c2 in this function to
    determine the first word. Returns boolean explaining if
    those concept pairs are compatible (do not overlap) """

    concepts_overlap_flag = False
    for c1_location in c1_dict_locations:
        for c2_location in c2_dict_locations:

            split_location_c1 = re.split(":", c1_location)
            split_location_c2 = re.split(":", c2_location)

            c1_start = int(split_location_c1[0])
            c1_end = int(split_location_c1[0]) + int(split_location_c1[1])
            c2_start = int(split_location_c2[0])
            c2_end = int(split_location_c2[0]) + int(split_location_c2[1])

            if c1_start > c2_start:
                c1_start, c2_start = c2_start, c1_start
                c1_end, c2_end = c2_end, c1_end

            if c2_start < c1_end:
                concepts_overlap_flag = True

    return concepts_overlap_flag


@task()
def import_concepts():
    """
    START (functions are above)

    This is where the number of documents to prepopulate the relation app starts
    from. For each document we make concept dictionaries from each pubtator containing
    information to 1) check for overlapping words 2) filter out bad documents
    3) Make sure there are at least two concept_types so that pairs can be formed
    4) Do not use concepts without a unique identification (UID)
    5) Use the longest word as the representative "text" if there are multiple texts
    per UID
    """

    from collections import Counter

    for document in Document.objects.all():

        df = document.as_pubtator_annotation_df()
        # We need annotations with at least a UID and Source
        df.dropna(subset=('uid', 'source'), how='any', inplace=True)
        #df = df[df['source'].isin(['MESH', 'MEDIC', 'NCBI Gene'])]

        # (TODO) Inspect for , in IDs and duplicate rows

        # Remove duplicate (uid, source, [not text]) and drop location + offset
        df.drop(['location', 'offset'], axis=1, inplace=True)
        df.drop_duplicates(['uid', 'source', 'text'], inplace=True)
        uid_source_uniq_df = df.drop_duplicates(['uid', 'source'])

        # Put these into the DB
        for index, row in df.iterrows():

            # (TODO) Maybe make these bulk creates
            c, created = Concept.objects.get_or_create(id=row['uid'])
            ct, created = ConceptText.objects.get_or_create(concept_id=row['uid'], text=row['text'])
            cdr, created = ConceptDocumentRelationship.objects.get_or_create(concept_text=ct, document=document, stype= stype_hash.get(row['source'], None) )


@task()
def compute_relationships():

    """This method takes a document and a relation pair list and makes the
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
                    document_id = document_id,
                    relation_type = relation_type,
                    concept_1_id = uid_a,
                    concept_2_id = uid_b )
