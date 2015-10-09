from django.conf import settings

from .models import Paper, Sentence, Relation, Concept  # JENNIFER TODO, import all models to be populated
from ..document.models import Document
import os
import subprocess
import re
import sys
from collections import defaultdict


def create_relation_concepts():
    """
    # ::TODO::, must restart Python shell? <- look into this
    # makes concept lists for only documents that have at least two dictionaries that contain concepts
    """
    queryset_documents = Document.objects.all()[:50]  # This will be changed when we improve the document dashboard

    for document in queryset_documents:
        chemical_dict, gene_dict, disease_dict, pubtator_bioc = document.make_concept_lists()
        total_full_dicts = 0
        if chemical_dict:
            total_full_dicts += 1
        if gene_dict:
            total_full_dicts += 1
        if disease_dict:
            total_full_dicts += 1
        if total_full_dicts > 1:
            document.make_cgd_concepts(chemical_dict, gene_dict, disease_dict)
            print "making concept dictionary for one document"
            concept_dict_list = [chemical_dict, gene_dict, disease_dict]
            document.add_relation_pairs_to_database(concept_dict_list)


create_relation_concepts()


if __name__ == "__main__":
    main()
