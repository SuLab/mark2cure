from django.conf import settings

from .models import Paper, Sentence, Relation, Concept  # JENNIFER TODO, import all models to be populated
from ..document.models import Document, Pubtator
import os
import subprocess
import re
import sys
from collections import defaultdict
import itertools

from mark2cure.common.bioc import BioCReader
from mark2cure.common.formatter import bioc_writer, bioc_as_json


"""
# makes concept lists for only documents that have at least two dictionaries
that contain concepts. Must be done server side so that it can be
determined if there are non-overlapping concepts
(at least one pair, or "relation pair" per document).

"""

def check_for_overlaps(c1_dict_location, c2_dict_location):
    """ Takes the location from the concepts and checks to see
    if they overlap. Swaps c1 and c2 in this function to
    determine the first word. Returns boolean explaining if
    those concept pairs are compatible (do not overlap) """

    split_location_c1 = re.split(":", c1_dict_location)
    split_location_c2 = re.split(":", c2_dict_location)

    c1_start = int(split_location_c1[0])
    c1_end = int(split_location_c1[0]) + int(split_location_c1[1])
    c2_start = int(split_location_c2[0])
    c2_end = int(split_location_c2[0]) + int(split_location_c2[1])

    if c1_start > c2_start:
        c1_start, c2_start = c2_start, c1_start
        c1_end, c2_end = c2_end, c1_end

    concepts_overlap = False
    if c2_start < c1_end:
        concepts_overlap = True

    return concepts_overlap


def return_pubtator_dict(pubtator, pub_type):
    """ Input one pubtator, and it's type and get one dictionary. TODO: later we might
    want to return different dictionaries."""
    if pubtator.valid():
        pub_dict = {}
        r = BioCReader(source=pubtator.content)
        r.read()
        for doc_idx, document in enumerate(r.collection.documents):
            for passage_idx, passage in enumerate(document.passages):
                for annotation in r.collection.documents[doc_idx].passages[passage_idx].annotations:
                    try:
                        text = annotation.text
                        uid = annotation.infons[pub_type]
                        location = str(annotation.locations[0])
                    except:
                        continue

                    if pub_type == "NCBI Gene":
                        stype = "g"
                    elif pub_type == "MEDIC":
                        stype = "d"
                    else:
                        stype = "c"

                    if uid != None:
                        pub_dict[uid] = {
                            'text': text,
                            'stype': stype,
                            'location': location,
                            'uid': uid
                            }

        return pub_dict


def make_concept_dicts_from_pubtators(document):
    """
    input a document and return all three special pubtators.
    """
    pub_gene = Pubtator.objects.get(document=document, kind="GNormPlus")
    pub_disease = Pubtator.objects.get(document=document, kind="DNorm")
    pub_chemical = Pubtator.objects.get(document=document, kind="tmChem")

    # TODO, if we add species, etc. this might change
    gene_dict =  return_pubtator_dict(pub_gene, "NCBI Gene")
    disease_dict = return_pubtator_dict(pub_disease, "MEDIC")
    chemical_dict = return_pubtator_dict(pub_chemical, "MESH")

    return gene_dict, disease_dict, chemical_dict


queryset_documents = Document.objects.all()[0:50]
for document in queryset_documents:

    gene_dict, disease_dict, chemical_dict = make_concept_dicts_from_pubtators(document)
    concept_dict_list = [gene_dict, disease_dict, chemical_dict]

    relation_pair_list = []
    for tuple_item in list(itertools.combinations(concept_dict_list, 2)):
        concept1_dict = tuple_item[0]
        concept2_dict = tuple_item[1]
        if concept1_dict and concept2_dict:
            for c1 in concept1_dict:
                for c2 in concept2_dict:

                    concepts_overlap = check_for_overlaps(concept1_dict[c1]['location'], concept2_dict[c2]['location'])

                    if concepts_overlap != True:
                        relation_pair_list.append([ concept1_dict[c1], concept2_dict[c2] ] )

    print relation_pair_list
    if relation_pair_list != []:
        document.add_concepts_to_database(relation_pair_list)
        document.add_relation_pairs_to_database(relation_pair_list)


if __name__ == "__main__":
    main()
