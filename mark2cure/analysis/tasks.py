'''
    Pulls data from xml and creates an array for each user consisting of PMID,
    type, and annotation. Uses NLTK scoring metrics tools to determine
    precision, recall, and f-score. By including PMID in the hash, this version
    allows for examining user to user comparisons across multiple documents in the
    group. Averages by User in one shot, instead of an average of averages.
    Uses userid instead of user_name. Treats one of the users as the test set, the
    other user as the gold standard for each pairing.
'''

from django.contrib.auth.models import User
from django.conf import settings

from ..common.formatter import clean_df
from ..common.models import Group
from ..document.models import Document
from .models import Report
from . import synonyms_dict

from nltk.metrics import scores as nltk_scoring
# from ..common import celery_app as app

import pandas as pd
import networkx as nx

import itertools


def hashed_er_annotations_df(group_pk, compare_type=True):
    """Generate a Entity Recognition DataFrame with additional hash column
    """
    group = Group.objects.get(pk=group_pk)
    org_er_df = Document.objects.ner_df(document_pks=group.get_document_pks(), include_pubtator=False)
    er_df = clean_df(org_er_df)

    if compare_type:
        er_df['hash'] = er_df.document_pk.apply(str) + '_' + er_df.ann_type_idx.apply(str) + '_' + er_df.section_offset.apply(str) + '_' + er_df.length.apply(str)
    else:
        er_df['hash'] = er_df.document_pk.apply(str) + '_' + er_df.section_offset.apply(str) + '_' + er_df.length.apply(str)
    return er_df


def compute_pairwise(hashed_er_anns_df):
    """
        Returns pairwise comparision between users (uesr_a & user_b)
        that have completed similar documents
    """
    # Make user_pks unique
    userset = set(hashed_er_anns_df.user_id)

    inter_annotator_arr = []
    # For each unique user comparision, compute
    for user_a, user_b in itertools.combinations(userset, 2):
        # The list of document_pks that each user had completed
        user_a_set = set(hashed_er_anns_df[hashed_er_anns_df['user_id'] == user_a].document_pk)
        user_b_set = set(hashed_er_anns_df[hashed_er_anns_df['user_id'] == user_b].document_pk)

        # Only compare documents both users have completed
        pmid_set = user_a_set.intersection(user_b_set)

        # If user_a and user_b have completed shared PMID, compute comparisions
        if len(pmid_set) != 0:
            pmid_df = hashed_er_anns_df[hashed_er_anns_df['document_pk'].isin(pmid_set)]
            ref_set = set(pmid_df[pmid_df['user_id'] == user_a].hash)
            test_set = set(pmid_df[pmid_df['user_id'] == user_b].hash)

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
    """
        Merging User1 and User2 columns for the pairings since combi ensures that
        that users are paired with each other only once (no reverse pairing)
    """
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

    temp_df = pd.DataFrame(all_users_arr, columns=('user_id', 'pairings', 'total_f'))

    # Obtaining average f-score from user-merged data.
    avg_f_arr = []
    for group_idx, group in temp_df.groupby('user_id'):
        pairing_counts = group['pairings'].sum()
        avg_f_arr.append((
            group_idx,
            pairing_counts,
            group['total_f'].sum() / pairing_counts
        ))

    avg_user_f = pd.DataFrame(avg_f_arr, columns=('user_id', 'pairings', 'f-score'))
    avg_user_f.sort('f-score', ascending=False, inplace=True)
    return avg_user_f


# @app.task(bind=True, ignore_result=True,
#           max_retries=0, soft_time_limit=600,
#           acks_late=True, track_started=True,
#           expires=3600)
def generate_reports(self, group_pk):
    args = locals()

    group = Group.objects.get(pk=group_pk)
    hash_table_df = hashed_er_annotations_df(group.pk)

    inter_annotator_df = compute_pairwise(hash_table_df)
    Report.objects.create(
        group=group, report_type=Report.PAIRWISE,
        dataframe=inter_annotator_df, args=args)

    avg_user_f = merge_pairwise_comparisons(inter_annotator_df)
    Report.objects.create(group=group, report_type=Report.AVERAGE,
            dataframe=avg_user_f, args=args)


# @app.task(bind=True, ignore_result=True,
#           max_retries=0,
#           acks_late=True, track_started=True,
#           expires=3600)
def group_analysis(self):
    for group in Group.objects.all():
        if group.percentage_complete() < 100:
            generate_reports.apply_async(
                args=[group.pk],
                queue='mark2cure_tasks')

    if not self.request.called_directly:
        return True


def hashed_annotations_graph_process(group_pk, min_thresh=settings.ENTITY_RECOGNITION_K):
    df = hashed_er_annotations_df(group_pk)

    # (TODO) This can be wayyy faster and better
    # Add a username column
    user_lookup = User.objects.all().values('pk', 'username')
    res = {}
    for u_dict in user_lookup:
        res[u_dict['pk']] = u_dict['username']

    df['username'] = df['user_id'].map(lambda user: res[user] if user > 0 else 'pubtator').apply(str)

    # Capitalize all annotation text
    df['text'] = df['text'].map(lambda x: x.upper())
    # Hard coded synonym cleaner
    df['clean_text'] = df['text'].map(lambda text_str: str(synonyms_dict.get(text_str, text_str))).apply(str)

    # Add field to deterine if hash meets minimum count
    hash_count_series = df['hash'].value_counts()
    df['hash_count'] = df['hash'].map(lambda hash_str: hash_count_series[hash_str])

    # Count the unique usage of that text string
    text_count_series = df['text'].value_counts()
    df['text_count'] = df['text'].map(lambda text_str: text_count_series[text_str])
    # clean_text_count_series = df['clean_text'].value_counts()
    # df['clean_text_count'] = df['clean_text'].map(lambda text_str: clean_text_count_series[text_str])

    # User Annotation count per PMID
    df['user_pmid_hash'] = df.document_pk.apply(str) + '_' + df.username.apply(str)
    user_pmid_hash_count_series = df['user_pmid_hash'].value_counts()
    df['user_pmid_count'] = df['user_pmid_hash'].map(lambda hash_str: user_pmid_hash_count_series[hash_str])

    return df[df['hash_count'] >= min_thresh]


def generate_network(group_pk: int, parallel=False, spring_force=10, include_degree=False): # noqa
    """
        1) Generate the DF needed to compute the Graph
        2) Compute the graph (aggregate, compare text / pmid)
        3) Compute graph metadata and attributes
    Args:
        group_pk (int): Use the group for selecting the ner documents to include

    Returns:

    """

    # Generate the required base DataFrame from raw Annotations
    df = hashed_annotations_graph_process(group_pk)

    # Numpy Arr of Unique Annotations via sanitized text
    nd_arr = df.clean_text.unique()
    # Unique Node labels (not using text as Identifier)
    names = ['n' + str(x + 1) for x in range(len(nd_arr))]
    # ID >> Cleantext lookup dictionary
    nodes = pd.Series(names, index=nd_arr).to_dict()

    # Start the Network off by adding all the unique Nodes (text annotations)
    G = nx.MultiGraph()
    G.add_nodes_from(names)

    anns = []
    for (doc, text, user), g_df in df.groupby(['document_pk', 'clean_text', 'username']):
        anns.append({
            'node': nodes[text],
            'doc': int(doc),
            'text': text,
            'user': user
        })
    new_df = pd.DataFrame(anns)

    if parallel:
        edge_idx = 1
        for (doc, user), g_df in new_df.groupby(['doc', 'user']):
            for node1, node2 in itertools.combinations(list(g_df.node), 2):
                G.add_edge(node1, node2, key=edge_idx, attributes={'doc': int(doc), 'user': user})
                edge_idx += 1

    else:
        edge_idx = 1
        for doc_id, group_df in new_df[['doc', 'node']].groupby('doc'):
            node_weight = group_df['node'].value_counts().T.to_dict()

            for node1, node2 in itertools.combinations(node_weight.keys(), 2):
                G.add_edge(node1, node2, key=edge_idx, attributes={'doc': int(doc_id)})
                edge_idx += 1

    # Compute Node Values
    pos = nx.spring_layout(G, iterations=spring_force)

    # Santize the node colors
    # if these colors get changed, need to edit group_home.js
    type_to_color = {0: '#d1f3ff', 1: '#B1FFA8', 2: '#ffd1dc'}
    df['color'] = df['ann_type_idx'].map(type_to_color)
    df['nodes'] = df['clean_text'].map(nodes)
    colors = {}
    for node, g_df in df.groupby('nodes'):
        colors[node] = g_df['color'].value_counts().idxmax()
    nx.set_node_attributes(G, 'color', colors)

    # Calcuate position and size
    x_pos, y_pos = {}, {}
    for idx, val in pos.items():
        x_pos[idx], y_pos[idx] = val
    nx.set_node_attributes(G, 'x', x_pos)
    nx.set_node_attributes(G, 'y', y_pos)

    size_dict = new_df['node'].value_counts().to_dict()
    for key, value in size_dict.items():
        size_dict[key] = int(value)

    nx.set_node_attributes(G, 'size', size_dict)
    nx.set_node_attributes(G, 'label', dict(zip(nodes.values(), nodes.keys())))

    if include_degree:
        # Calculate centrality metrics
        degree = nx.degree_centrality(G)
        attributes = {}
        for key, value in degree.iteritems():
            attributes[key] = {}
            attributes[key]['degree'] = int(value)
        nx.set_node_attributes(G, 'attributes', attributes)

    return G



