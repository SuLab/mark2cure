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
from ..document.models import Section, Annotation
from ..task.entity_recognition.models import EntityRecognitionAnnotation
from .models import Report

from django.contrib.auth.models import User

from nltk.metrics import scores as nltk_scoring
from celery import task

import pandas as pd
import networkx as nx

import itertools


def hashed_annotations_df(group_pk, private_api=False,
        compare_type=True, compare_text=False):
    '''
        Users will always be represented in Analysis as their PKs as strings.
        This is because, on occasion, we will include users that are not in the
        DB and don't have an int to represent their PK

        This standard is a must as there are many downstream cases where this
        column is searched for users by their str(pk)
    '''

    hashed_annotations = []  # Creates empty list for the annotations
    ann_types = EntityRecognitionAnnotation.ANNOTATION_TYPE_CHOICE

    if private_api:
        group = Group.objects.get(pk=group_pk)
        document_pmids = group.get_documents().values_list('document_id', flat=True)

        # Fetch all the section pks and their text length
        sections = Section.objects.filter(
            document__document_id__in=document_pmids
        ).extra(select={'section_length': 'LENGTH(text)'}).values('pk', 'section_length')

        for doc_pmid in document_pmids:

            document_annotations = EntityRecognitionAnnotation.objects.annotations_for_document_pmid(doc_pmid)

            passage_offset = 0

            section_pks = list(set([ann.section_id for ann in document_annotations]))
            for section_pk in section_pks:

                section_annotations = filter(lambda ann: ann.section_id == section_pk, document_annotations)
                for ann in section_annotations:
                    hashed_annotations.append((
                        ann.document_id,
                        str(ann.user_id),
                        ann_types.index(ann.type),

                        passage_offset + ann.start,

                        len(ann.text),
                        ann.text
                    ))

                section_results = filter(lambda section: section['pk'] == section_pk, sections)
                passage_offset += int(section_results[0]['section_length'])

    else:
        # Capture exported data
        req = group_users_bioc({}, group_pk, 'xml')
        reader = BioCReader(source=req.content)
        reader.read()

        # Read through BioC results and convert to a list of (user, uniq_ann_identifier, document_id)
        for document in reader.collection.documents:
            for passage in document.passages:
                for ann in passage.annotations:
                    ann_loc = ann.locations[0]

                    hashed_annotations.append((
                        document.id,
                        str(ann.infons.get('user')),
                        ann_types.index(ann.infons.get('type')),
                        ann_loc.offset,
                        ann_loc.length,
                        ann.text
                    ))

    df = pd.DataFrame(hashed_annotations, columns=('document_id', 'user', 'type', 'offset', 'length', 'text'))

    if compare_type:
        df['hash'] = df.document_id.apply(str) + '_' + df.type.apply(str) + '_' + df.offset.apply(str) + '_' + df.length.apply(str)
    else:
        df['hash'] = df.document_id.apply(str) + '_' + df.offset.apply(str) + '_' + df.length.apply(str)

    return df


def compute_pairwise(hashed_annotations_df):
    '''
        Returns pairwise comparision between users (uesr_a & user_b)
        that have completed similar documents
    '''
    # Make user_pks unique
    userset = set(hashed_annotations_df.user)

    inter_annotator_arr = []
    # For each unique user comparision, compute
    for user_a, user_b in itertools.combinations(userset, 2):
        # The list of document_ids that each user had completed
        user_a_set = set(hashed_annotations_df[hashed_annotations_df['user'] == user_a].document_id)
        user_b_set = set(hashed_annotations_df[hashed_annotations_df['user'] == user_b].document_id)

        # Only compare documents both users have completed
        pmid_set = user_a_set.intersection(user_b_set)

        # If user_a and user_b have completed shared PMID, compute comparisions
        if len(pmid_set) != 0:
            pmid_df = hashed_annotations_df[hashed_annotations_df['document_id'].isin(pmid_set)]
            ref_set = set(pmid_df[pmid_df['user'] == user_a].hash)
            test_set = set(pmid_df[pmid_df['user'] == user_b].hash)

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
    hash_table_df = hashed_annotations_df(group_pk, private_api=private_api)

    inter_annotator_df = compute_pairwise(hash_table_df)
    Report.objects.create(
        group=group, report_type=Report.PAIRWISE,
        dataframe=inter_annotator_df, args=args)

    avg_user_f = merge_pairwise_comparisons(inter_annotator_df)
    Report.objects.create(group=group, report_type=Report.AVERAGE,
            dataframe=avg_user_f, args=args)


def hashed_annotations_graph_process(group_pk, min_thresh=15):
    df = hashed_annotations_df(group_pk, private_api=True)

    # (TODO) This can be wayyy faster and better
    # Add a username column
    user_lookup = User.objects.all().values('pk', 'username')
    res = {}
    for u_dict in user_lookup:
        res[u_dict['pk']] = u_dict['username']
    df['username'] = df['user'].map(lambda user: res[int(user)]).apply(str)

    # Capitalize all annotation text
    df['text'] = df['text'].map(lambda x: x.upper())
    # Hard coded synonym cleaner
    synonyms_dict = pd.read_csv('scripts/synonym_dictionary.txt', sep='\t', names=['dirty', 'clean'], index_col='dirty').to_dict()['clean']
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
    df['user_pmid_hash'] = df.document_id.apply(str) + '_' + df.user.apply(str)
    user_pmid_hash_count_series = df['user_pmid_hash'].value_counts()
    df['user_pmid_count'] = df['user_pmid_hash'].map(lambda hash_str: user_pmid_hash_count_series[hash_str])

    return df[df['hash_count'] >= min_thresh]


def generate_network(group_pk, parallel=False, spring_force=10, include_degree=False):
    '''
        1) Generate the DF needed to compute the Graph
        2) Compute the graph (aggregate, compare text / pmid)
        3) Compute graph metadata and attributes
    '''

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
    for (doc, text, user), g_df in df.groupby(['document_id', 'clean_text', 'username']):
        anns.append({
            'node': nodes[text],
            'doc': doc,
            'text': text,
            'user': user
        })
    new_df = pd.DataFrame(anns)

    if parallel:
        edge_idx = 1
        for (doc, user), g_df in new_df.groupby(['doc', 'user']):
            for node1, node2 in itertools.combinations(list(g_df.node), 2):
                G.add_edge(node1, node2, key=edge_idx, attributes={'doc': doc, 'user': user})
                edge_idx += 1

    else:
        edge_idx = 1
        for doc_id, group_df in new_df[['doc', 'node']].groupby('doc'):
            node_weight = group_df['node'].value_counts().T.to_dict()

            for node1, node2 in itertools.combinations(node_weight.keys(), 2):
                G.add_edge(node1, node2, key=edge_idx, attributes={'doc': doc_id})
                edge_idx += 1

    # Compute Node Values
    pos = nx.spring_layout(G, iterations=spring_force)

    # Santize the node colors
    type_to_color = {0: '#d1f3ff', 1: '#B1FFA8', 2: '#ffd1dc'}
    df['color'] = df['type'].map(type_to_color)
    df['nodes'] = df['clean_text'].map(nodes)
    colors = {}
    for node, g_df in df.groupby('nodes'):
        colors[node] = g_df['color'].value_counts().idxmax()
    nx.set_node_attributes(G, 'color', colors)

    # Calcuate position and size
    x_pos, y_pos = {}, {}
    for idx, val in pos.iteritems():
        x_pos[idx], y_pos[idx] = val
    nx.set_node_attributes(G, 'x', x_pos)
    nx.set_node_attributes(G, 'y', y_pos)

    nx.set_node_attributes(G, 'size', new_df['node'].value_counts().to_dict())
    nx.set_node_attributes(G, 'label', dict(zip(nodes.values(), nodes.keys())))

    if include_degree:
        # Calculate centrality metrics
        degree = nx.degree_centrality(G)
        attributes = {}
        for key, value in degree.iteritems():
            attributes[key] = {}
            attributes[key]['degree'] = value
        nx.set_node_attributes(G, 'attributes', attributes)

    return G


@task()
def group_analysis():
    for group in Group.objects.all():
        if group.percentage_complete() < 100:
            generate_reports(group.pk, private_api=True)
