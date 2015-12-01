from ...common.bioc import BioCReader
from ...document.models import Document, Pubtator

import re
import itertools

"""
# makes concept lists for only documents that have at least two dictionaries
that contain concepts. Must be done server side so that it can be
determined if there are non-overlapping concepts
(at least one pair, or "relation pair" per document).
"""

bad_document_flag = False


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
    global bad_document_flag
    if pubtator.valid():
        pub_dict = {}
        r = BioCReader(source=pubtator.content)
        r.read()
        for doc_idx, document in enumerate(r.collection.documents):
            for passage_idx, passage in enumerate(document.passages):

                if bad_document_flag is True:
                    # flag already raised, so don't need to check for the following
                    pass
                elif "%" in passage.text or "<" in passage.text or ">" in passage.text:
                    bad_document_flag = True

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
queryset_documents = Document.objects.all()[:500]

for document in queryset_documents:
    concept_dict_list = make_concept_dicts_from_pubtators(document)
    """list of docs that we know have bad Pubtators (they'll produce bad
    highlights). They contain characters that throw of the Pubtator's results
    such as >, < or %. Causes unpredictable results.
    Ex: document pk 1801, Sometimes percent signs
    are stored as percents in the M2C DB, but sometimes they come out as a "y"
    like in 1801.
    """
    bad_document_list = ["1801"]

    if bad_document_flag is True or str(document.pk) in bad_document_list:
        # reset flag and continue through documents
        bad_document_flag = False
        continue
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
        document.add_concepts_to_database(relation_pair_list)
        document.add_relation_pairs_to_database(relation_pair_list)


if __name__ == "__main__":
    main()
