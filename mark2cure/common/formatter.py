from ..common.bioc import BioCWriter, BioCCollection, BioCAnnotation, BioCLocation
from ..document.models import Document

import xmltodict
import itertools
import datetime
import json
import nltk


def pad_split(text):
    text = text.replace("\\(", " ( ")
    text = text.replace("\\)", " ) ")
    text = text.replace("\\.", " . ")
    text = text.replace("\\,", " , ")
    text = text.replace("\\%", " % ")
    text = text.replace("\\#", " # ")
    text = text.replace("\\&", " & ")
    text = text.replace("\\+", " + ")
    text = text.replace("\\=", " = ")
    text = text.replace("\\[", " [ ")
    text = text.replace("\\]", " ] ")
    text = text.replace("\\;", " ; ")
    text = text.replace("\\/", " / ")
    text = text.replace("/", " / ")
    text = text.replace("\\\"", " \" ")
    text = text.replace("  ", " ")
    text = text.replace("  ", " ")
    return nltk.word_tokenize(text.encode('utf-8'))


def bioc_writer(request):
    writer = BioCWriter()
    writer.collection = BioCCollection()
    writer.collection.date = datetime.date.today().strftime("%Y%m%d")
    if request:
        writer.collection.source = 'Mark2Cure API: {relative_url}'.format(
            relative_url=request.META.get('PATH_INFO', ''))
    else:
        writer.collection.source = 'Mark2Cure Internal'
    return writer


def bioc_as_json(writer):
    o = xmltodict.parse(writer.__str__())

    try:
        json_dict = json.loads(json.dumps(o))
        for passage_index, passage in enumerate(json_dict.get('collection').get('document').get('passage')):

            anns = passage.get('annotation')
            if type(anns) != list:
                json_dict['collection']['document']['passage'][passage_index]['annotation'] = [anns]

        return json.dumps(json_dict)

    except:
        return json.dumps(o)


# R & S are tuple of (start position, stop position)
def are_separate(r, s):
    return r[1] < s[0] or s[1] < r[0]


def are_overlapping(r, s):
    return not(r[1] < s[0] or s[1] < r[0])


def clean_passage_df_overlaps(df):
    df.reset_index(inplace=True)

    res = []
    for index, row in df.iterrows():
        loc_splt_str = [int(x) for x in row['location'].split(':')]
        loc_span = (loc_splt_str[0], loc_splt_str[0] + loc_splt_str[1])
        res.append((loc_span, index))

    for x, y in itertools.combinations(res, 2):
        span_a, span_a_row = x
        span_b, span_b_row = y
        if are_overlapping(span_a, span_b):
            df.drop(span_b_row, inplace=True)

    return df


def clean_df(df):
    ann_types_arr = ['Chemical', 'Gene', 'Disease']
    ann_types_arr.extend(Document.APPROVED_TYPES)
    # (TODO) Consolidate into a single function for reuse with task.relation.task.import_concepts
    df.dropna(subset=('uid', 'source'), how='any', inplace=True)
    df = df[df['ann_type'].isin(ann_types_arr)]

    # (TODO) Inspect for , in IDs and duplicate rows
    # remove unnecessary prefixes from uids
    df.loc[:, "uid"] = df.loc[:, "uid"].map(lambda v: v[5:] if v.startswith("MESH:") else v)
    df.loc[:, "uid"] = df.loc[:, "uid"].map(lambda v: v[5:] if v.startswith("OMIM:") else v)
    df.loc[:, "uid"] = df.loc[:, "uid"].map(lambda v: v[6:] if v.startswith("CHEBI:") else v)

    # (TODO) Is there an ordering to the UIDs?
    # (NOTES) After a short inspection, I didn't see an obvious order. -Max 3/2/2016
    df = df[~df.uid.str.contains(",")]
    df = df[~df.uid.str.contains("\|")]

    return df


def apply_annotations(writer, df):
    '''
        This takes a BioCWriter for 1 document and a Pandas Dataframe of annotations
        to apply to the BioCWriter
    '''

    # If nothing in DF, we can safely return the non-modified writer
    if df.shape[0] == 0:
        return writer

    # (TODO) Cases in which a user hasn't annotated something in the title...
    section_ids = list(df['section_id'].unique())
    if not len(section_ids) >= 1:
        raise ValueError('Incorrect number of document sections.')

    bioc_doc = writer.collection.documents[0]
    approved_types = ['disease', 'gene_protein', 'chemical',    'gene', 'species']
    # ann_types_arr = ['Chemical', 'Gene', 'Disease']
    # ann_types_arr.extend(Document.APPROVED_TYPES)

    # Make all the offsets scoped the the entire document (like Pubtator)
    df.ix[df['offset_relative'], 'start_position'] = df['section_offset'] + df['start_position']
    df.ix[df['offset_relative'], 'offset_relative'] = False

    # Not required, but easier to view this way
    df.sort('start_position', inplace=True)

    # if row['ann_type'] in approved_types:

    i = 0
    for offset_position, group_df in df.groupby(['section_offset']):
        # This is another approach to offset grouping:
        # offset = int(bioc_passage.offset)
        # for idx, row in df[df['offset'] == offset].iterrows():

        bioc_passage = bioc_doc.passages[i]
        i = i + 1

        # Should already be empty if from doc.as_writer(), but perhaps we add an
        # append method in the future
        bioc_passage.clear_annotations()

        # (TODO) This should determined off the user_id field
        # if pubtator:
        #     group_df = clean_passage_df_overlaps(group_df)

        for row_idx, row in group_df.iterrows():
            annotation = BioCAnnotation()

            annotation.id = str(row_idx)

            # if pubtator:
            #    annotation.put_infon('user', 'pubtator')

            annotation.put_infon('uid', row['uid'])
            annotation.put_infon('source', row['source'])
            annotation.put_infon('user_id', str(row['user_id']))
            annotation.put_infon('type', row['ann_type'])
            annotation.put_infon('type_id', str(approved_types.index(row['ann_type'])))

            location = BioCLocation()
            location.offset = str(row['start_position'])
            location.length = str(row['length'])
            annotation.add_location(location)

            annotation.text = row['text']

            bioc_passage.add_annotation(annotation)

    return writer

