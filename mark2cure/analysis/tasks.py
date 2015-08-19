'''
    Pulls data from xml and creates an array for each user consisting of PMID,
    type, and annotation. Uses NLTK scoring metrics tools to determine
    precision, recall, and f-score. By including PMID in the hash, this version
    allows for examining user to user comparisons across multiple documents in the
    group. Averages by User in one shot, instead of an average of averages.
    Uses userid instead of user_name. Treats one of the users as the test set, the
    other user as the gold standard for each pairing.
'''

from ..common.bioc import BioCReader
from ..common.models import Group
from ..api.views import group_users_bioc
from ..document.models import Annotation
from .models import Report

from nltk.metrics import scores as nltk_scoring
from celery import task

import pandas as pd
import itertools


def hashed_annotations_df(group_pk, private_api=False,
        compare_type=True, compare_text=False):

    if private_api:
        pass

    else:
        # Capture exported data
        req = group_users_bioc({}, group_pk, 'xml')
        reader = BioCReader(source=req.content)
        reader.read()

        hashed_annotations = [] # Creates empty list for the annotations
        ann_types = Annotation.ANNOTATION_TYPE_CHOICE

        # Read through BioC results and convert to a list of (user, uniq_ann_identifier, document_id)
        for document in reader.collection.documents:
            for passage in document.passages:
                for ann in passage.annotations:
                    ann_loc = ann.locations[0]

                    hashed_annotations.append((
                        document.id,
                        ann.infons.get('user'),
                        ann_types.index( ann.infons.get('type') ),
                        ann_loc.offset,
                        ann_loc.length,
                        ann.text
                    ))

    df = pd.DataFrame(hashed_annotations, columns=('document_id', 'user', 'type', 'offset', 'length', 'text'))

    if compare_type:
        df['hash'] = df.document_id.apply(str) +'_'+ df.type.apply(str) +'_'+ df.offset.apply(str) +'_'+ df.length.apply(str)
    else:
        df['hash'] = df.document_id.apply(str) +'_'+                          df.offset.apply(str) +'_'+ df.length.apply(str)

    return df


def compute_pairwise(hashed_annotations_df):
    '''
        Returns pairwise comparision between users (uesr_a & user_b)
        that have completed similar documents
    '''
    # Make user_pks unique
    userset = set(hashed_annotations_df.user)

    # (TODO) Why does this need to be sorted by User
    # an_srtdby_user = sorted(hashed_annotations)
    # inter_annotator_df.sort('document_id', inplace=True)

    inter_annotator_arr = []
    # For each unique user comparision, compute
    for user_a, user_b in itertools.combinations(userset, 2): # (TODO) Put counter in here
        # The list of document_ids that each user had completed
        user_a_set = set(hashed_annotations_df[hashed_annotations_df['user'] == user_a].document_id)
        user_b_set = set(hashed_annotations_df[hashed_annotations_df['user'] == user_b].document_id)

        # Only compare documents both users have completed
        pmid_set = user_a_set.intersection(user_b_set)

        # If user_a and user_b have completed shared PMID, compute comparisions
        if len(pmid_set) != 0:
            pmid_df = hashed_annotations_df[hashed_annotations_df['document_id'].isin(pmid_set)]
            ref_set  = set(pmid_df[ pmid_df['user'] == user_a ].hash)
            test_set = set(pmid_df[ pmid_df['user'] == user_b ].hash)

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

    return pd.DataFrame(inter_annotator_arr, columns=('user_a', 'user_b', 'docs_compared', 'precision', 'recall', 'f-score'))


def merge_pairwise_comparisons(inter_annotator_df):
    '''
        Merging User1 and User2 columns for the pairings since combi ensures that
        that users are paired with each other only once (no reverse pairing)
    '''
    # Sort rows by best F-Score at the top
    inter_annotator_df.sort('f-score', ascending=False, inplace=True)

    all_users_arr = []
    for group_idx, group in inter_annotator_df.groupby('user_a'):
        all_users_arr.append((
            group_idx,
            group.shape[0],
            group['f-score'].sum()
        ))

    for group_idx, group in inter_annotator_df.groupby('user_b'):
        all_users_arr.append((
            group_idx,
            group.shape[0],
            group['f-score'].sum()
        ))

    temp_df = pd.DataFrame(all_users_arr, columns=('user', 'pairings', 'total_f'))

    # Obtaining average f-score from user-merged data.
    # print 'Obtaining average f-scores'
    avg_f_arr = []
    for group_idx, group in temp_df.groupby('user'):
        pairing_counts = group['pairings'].sum()
        avg_f_arr.append((
            group_idx,
            pairing_counts,
            group['total_f'].sum() / pairing_counts
        ))

    avg_user_f = pd.DataFrame(avg_f_arr, columns=('user', 'pairings', 'f-score'))
    avg_user_f.sort('f-score', ascending=False, inplace=True)
    return avg_user_f


def generate_reports(group_pk, private_api=False,
        compare_type=True, compare_text=False):
    args = locals()

    group = Group.objects.get(pk=group_pk)
    hash_table_df = hashed_annotations_df(group_pk)

    inter_annotator_df = compute_pairwise(hash_table_df)
    Report.objects.create(
            group=group, report_type=1,
            dataframe=inter_annotator_df, args=args)

    avg_user_f = merge_pairwise_comparisons(inter_annotator_df)
    Report.objects.create(group=group, report_type=1,
            dataframe=avg_user_f, args=args)


@task()
def group_analysis():
    for group_pk in Group.objects.values_list('pk', flat=True):
        generate_reports(group_pk)
