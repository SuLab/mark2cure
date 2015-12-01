from ...common.bioc import BioCReader
from ...document.models import Document, Pubtator
from .models import Concept, Relation, Answer

import re
import itertools
from celery import task

"""
# makes concept lists for only documents that have at least two dictionaries
that contain concepts. Must be done server side so that it can be
determined if there are non-overlapping concepts
(at least one pair, or "relation pair" per document).
"""



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


def return_pubtator_dict(pubtator, pub_key, stype):
    """ Input one pubtator (GNormPlus, DNorm, tmChem) and the pub type key (NCBI
    Gene, MEDIC, or MESH) for the infons. And the stype and return one
    dictionary. TODO: later we might want to return different dictionaries so
    we would use different 'pub keys' and stypes."""

    if pubtator.valid():
        pub_dict = {}
        r = BioCReader(source=pubtator.content)
        r.read()
        for doc_idx, document in enumerate(r.collection.documents):
            for passage_idx, passage in enumerate(document.passages):

                if "%" in passage.text or "<" in passage.text or ">" in passage.text:
                    raise ValueError('Failed to parse b/c of bad characters')

                for annotation in r.collection.documents[doc_idx].passages[passage_idx].annotations:
                    try:
                        text = annotation.text
                        uid = annotation.infons[pub_key]
                        location = str(annotation.locations[0])
                    except:
                        continue

                    if uid is not None:
                        if uid in pub_dict:
                            if location not in pub_dict[uid]['location']:
                                pub_dict[uid]['location'].append(location)
                            if text not in pub_dict[uid]['text']:
                                pub_dict[uid]['text'].append(text)
                        else:
                            pub_dict[uid] = {
                                'text': [text],
                                'stype': stype,
                                'location': [location],
                                'uid': uid
                            }
        return pub_dict


def make_concept_dicts_from_pubtators(document):
    """Input a document and return information from all three Pubtators (if
    there are three). Try/except because sometimes there are no pubtators.
    """
    available_pub_list = []
    try:
        pub_gene = Pubtator.objects.get(document=document, kind="GNormPlus")
        gene_dict = return_pubtator_dict(pub_gene, "NCBI Gene", 'g')
        available_pub_list.append(gene_dict)
    except:
        print document.document_id, "CHECK THE PUBTATOR gnorm"
        pass

    try:
        pub_disease = Pubtator.objects.get(document=document, kind="DNorm")
        disease_dict = return_pubtator_dict(pub_disease, "MEDIC", 'd')
        available_pub_list.append(disease_dict)
    except:
        print document.document_id, "CHECK THE PUBTATOR dnorm"
        pass

    try:
        pub_chemical = Pubtator.objects.get(document=document, kind="tmChem")
        chemical_dict = return_pubtator_dict(pub_chemical, "MESH", 'c')
        available_pub_list.append(chemical_dict)
    except:
        print document.document_id, "CHECK THE PUBTATOR, tmchem"
        pass

    # TODO, if we add species, etc. this might change
    pub_dict_list = []
    for dict_item in available_pub_list:
        if dict_item:
            pub_dict_list.append(dict_item)

    return pub_dict_list

def add_concepts_to_database(document, relation_pair_list):
    """This method takes a document and relation pair list and makes
    annotations/unique concepts for the document. If annotations already
    exist for this document, then no annotations will be made. Concepts are
    only created when a variety of criteria is met (see relation/tasks for
    a description of the criteria)
    """
    if not Concept.objects.filter(document=document).exists():
        unique_id_list = []  # avoid duplicates
        for pair in relation_pair_list:
            for concept in pair:
                if concept['uid'] in unique_id_list:
                    continue
                else:
                    unique_id_list.append(concept['uid'])
                    Concept.objects.create(document=document, uid=concept['uid'], stype=concept['stype'], text=max(concept['text'], key=len) )


def add_relation_pairs_to_database(document, relation_pair_list):
    """This method takes a document and a relation pair list and makes the
    appropriate "relation." The reason relations are stored instead of
    created on the fly, is because there was a significant amount of code
    performed to determine whether or not relations even exist.  We need
    those relation pairs to determine if the document is worth showing to
    a user.
    """
    if not Relation.objects.filter(document=document).exists():

        for pair in relation_pair_list:
            concept1 = pair[0]
            concept2 = pair[1]
            relation_type = concept1['stype'] + "_" + concept2['stype']

            # always have gene on the left and disease on the right for
            # sentence structure (this is important for future analysis)
            # Keeps concepts inside of "relations" in right orientation in db
            # TODO (IMPORTANT if adding more concepts)
            if relation_type == 'c_g' or relation_type == 'd_g' or relation_type == 'd_c':
                concept1, concept2 = concept2, concept1

            concept1_obj = Concept.objects.get(document=document, uid=concept1['uid'])
            concept2_obj = Concept.objects.get(document=document, uid=concept2['uid'])
            Relation.objects.create(document=document, relation=[ concept1['uid'], concept2['uid'] ], concept1_id=concept1_obj, concept2_id=concept2_obj )


@task()
def import_document_relationships():
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

    # this is how many documents to import via ***TASKS***
    queryset_documents = Document.objects.all()[:5000]

    for document in queryset_documents:
        concept_dict_list = make_concept_dicts_from_pubtators(document)
        # (TODO) This needs to be fixed, can't just have a bad doc list to avoid -Max
        """list of docs that we know have bad Pubtators (they'll produce bad
        highlights). They contain characters that throw of the Pubtator's results
        such as >, < or %. Causes unpredictable results.
        Ex: document pk 1801, Sometimes percent signs
        are stored as percents in the M2C DB, but sometimes they come out as a "y"
        like in 1801.
        """

        print '-'*30
        print concept_dict_list

        relation_pair_list = []
        if len(concept_dict_list) >= 2:
            for tuple_item in list(itertools.combinations(concept_dict_list, 2)):
                concept1_dict = tuple_item[0]
                concept2_dict = tuple_item[1]
                if concept1_dict and concept2_dict:
                    for c1 in concept1_dict:
                        for c2 in concept2_dict:

                            concepts_overlap_flag = check_for_overlaps(concept1_dict[c1]['location'], concept2_dict[c2]['location'])
                            """Maximum is the longest word out of all possible
                            "texts" for one unique ID (preference to show word and
                            not the acronym
                            """
                            if concepts_overlap_flag is not True and max(concept1_dict[c1]['text'], key=len) != max(concept2_dict[c2]['text'], key=len):
                                relation_pair_list.append([concept1_dict[c1], concept2_dict[c2]])

        """Add the concepts for each document, then add the relations for each
        document
        """
        if relation_pair_list != []:
            add_concepts_to_database(document, relation_pair_list)
            add_relation_pairs_to_database(document, relation_pair_list)

