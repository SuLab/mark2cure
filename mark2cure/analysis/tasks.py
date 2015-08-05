'''
    Pulls data from xml and creates an array for each user consisting of PMID,
    type, and annotation. Uses NLTK scoring metrics tools to determine
    precision, recall, and f-score. By including PMID in the hash, this version
    allows for examining user to user comparisons across multiple documents in the
    group. Averages by User in one shot, instead of an average of averages.
    Uses userid instead of user_name. Treats one of the users as the test set, the
    other user as the gold standard for each pairing.
'''

from operator import itemgetter
from ..common.bioc import BioCReader

from nltk.metrics import scores as nltk_scoring
import pandas as pd

import itertools
import requests
import numpy
import os


def inter_annotator_agreement(group_id):

    # Capture exported data
    req = requests.get('https://mark2cure.org/api/group/'+str(group_id)+'/user_annotations.xml')
    assert req.ok, "BioC Fetch Failed"
    reader = BioCReader(source=req.content)
    reader.read()

    HashedInfo = [] #Creates empty list for the annotations
    user_pks = []

    # Read through BioC results and convert to a list of (user, uniq_ann_identifier, document_id)
    for document in reader.collection.documents:
        # print "for PMID: ",document.id
        for passage in document.passages:
            for ann in passage.annotations:
                infons = ann.infons
                user_pks.append(str(ann.infons.get('user')))
                ann_loc = ann.locations[0]

                HashedInfo.append((
                    infons.get('user'),
                    '{0}_{1}_{2}_{3}'.format(document.id, infons.get('type'), ann_loc.offset, ann_loc.length),
                    document.id ))

    # Make user_pks unique
    userset = set(user_pks)

    print "generating user pairs and calculating precision, recall, and f"
    # Sorts the the list
    an_srtdby_user = sorted(HashedInfo)

    inter_annotator_arr = []
    # For each unique user comparision, compute
    for user_a, user_b in itertools.combinations(userset, 2): # (TODO) Put counter in here
        ref_set = set()  # Empties the set prior to each match
        test_set = set()  # Empties the set prior to each match
        user_a_set = set()  # Empties the set prior to each match
        user_b_set = set()  # Empties the set prior to each match
        pmid_set = set()  # Empties the set prior to each match

        for user_pk, uniq_hash, pmid in an_srtdby_user:
            if user_pk == user_a:
                user_a_set.add(pmid)

            if user_pk == user_b:
                user_b_set.add(pmid)

        # Only compare documents both users have completed
        pmid_set = user_a_set.intersection(user_b_set)
        for ann_pmid in pmid_set:
            for user_pk, uniq_hash, pmid in an_srtdby_user:

                if pmid == ann_pmid:
                    if user_pk == user_a:
                        ref_set.add(uniq_hash)

                    if user_pk == user_b:
                        test_set.add(uniq_hash)

        # If user_a and user_b have completed shared PMID, compute comparisions
        if len(pmid_set) != 0:
            # Compute the precision, recall and F-measure based on
            # the unique hashes
            inter_annotator_arr.append((
                user_a,
                user_b,
                len(pmid_set),
                nltk_scoring.precision(ref_set, test_set),
                nltk_scoring.recall(ref_set, test_set),
                nltk_scoring.f_measure(ref_set, test_set)
            ))

    inter_annotator_df = pd.DataFrame(inter_annotator_arr, columns=('user_a', 'user_b', 'docs_compared', 'precision', 'recall', 'f-score'))
    inter_annotator_df.to_csv('max_documentset_fscore.csv', index=False)

    all_users_arr = []
    ### Merging User1 and User2 columns for the pairings since combi ensures that
    ### that users are paired with each other only once (no reverse pairing)
    for groupByuser1, rows in itertools.groupby(inter_annotator_df.values, key=itemgetter(0)):
        grpd_a, counter = 0, 0
        for user_a, user_b, docs_compared, precision, recall, fscore in rows:
            grpd_a += fscore
            counter += 1

        # dumps total f-score and total counts for each user
        all_users_arr.append((
            groupByuser1,
            counter,
            grpd_a
        ))

    for groupByuser2, rows in itertools.groupby(inter_annotator_df.values, key=itemgetter(1)):
        grpd_b, counter = 0, 0
        for user_a, user_b, docs_compared, precision, recall, fscore in rows:
            grpd_b += fscore
            counter += 1

        # dumps total f-score and total counts for each user
        all_users_arr.append((
            groupByuser2,
            counter,
            grpd_a
        ))
    avg_user_f_sum = pd.DataFrame(all_users_arr, columns=('user', 'pairings', 'avg_f'))
    avg_user_f_sum.to_csv('max_avg_fscore_sum.csv', index=False)

    # Load the Avg user F-Score data and sort by F
    al_user_data = numpy.sort(avg_user_f_sum.values, axis=0)

    # Obtaining average f-score from user-merged data.
    print 'Obtaining average f-scores'
    avg_f_arr = []
    for groupByusers, rows in itertools.groupby(al_user_data, key=itemgetter(0)):
        grpd_counts, grpd_c, counter = 0, 0, 0

        for user_pk, pairings, avg_f in rows:
            grpd_c += avg_f
            grpd_counts += pairings
            counter += 1

        avg_f_arr.append((
            groupByusers,
            grpd_counts,
            grpd_c/grpd_counts
        ))

    avg_user_f = pd.DataFrame(avg_f_arr, columns=('user', 'pairings', 'avg_f'))
    avg_user_f.to_csv('max_avg_fscore.csv', index=False)
