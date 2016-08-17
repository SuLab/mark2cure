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


def is_pubtator_df(df):
    return df['user_id'].isnull().all()


def clean_df(df, overlap_protection=False):
    '''
        This attempts to santize our Annotation Dataframes that may originate
        from multiple sources (users, pubtator) so they're comparable
    '''
    if df.shape[1] is not 10:
        raise ValueError('Incorrect number of dataframe columns.')

    # We're previously DB Primary Keys
    df.reset_index(inplace=True)

    # Make all the offsets scoped the the entire document (like Pubtator)
    df.ix[df['offset_relative'], 'start_position'] = df['section_offset'] + df['start_position']
    df.ix[df['offset_relative'], 'offset_relative'] = False

    # Not required, but easier to view this way
    df.sort('start_position', inplace=True)

    # Absolutely require UID and Source
    df.dropna(subset=('uid', 'source'), how='any', inplace=True)

    # Remove unnecessary prefixes from uids if coming from external sources (via pubtator algos)
    df.loc[:, 'uid'] = df.loc[:, 'uid'].map(lambda v: v[5:] if v.startswith('MESH:') else v)
    df.loc[:, 'uid'] = df.loc[:, 'uid'].map(lambda v: v[5:] if v.startswith('OMIM:') else v)
    df.loc[:, 'uid'] = df.loc[:, 'uid'].map(lambda v: v[6:] if v.startswith('CHEBI:') else v)

    # (TODO) Inspect for , in IDs and duplicate rows
    # (TODO) Is there an ordering to the UIDs?

    # (NOTES) After a short inspection, I didn't see an obvious order. -Max 3/2/2016
    df = df[~df.uid.str.contains(",")]
    df = df[~df.uid.str.contains("\|")]

    # Only keep rows that are in our known annotation type lists
    df['ann_type'] = df['ann_type'].str.lower()
    ann_types_arr = ['chemical', 'gene', 'disease']
    ann_types_arr.extend(Document.APPROVED_TYPES)
    df = df[df['ann_type'].isin(ann_types_arr)]
    df['ann_type_id'] = 0

    # is_pubtator = is_pubtator_df(df)
    if overlap_protection:
        # Removes any annotations (rows) that have span overlap
        res = []
        for index, row in df.iterrows():
            loc_span = (row['start_position'], row['start_position'] + row['length'])
            res.append((loc_span, index))

        for x, y in itertools.combinations(res, 2):
            span_a, span_a_row = x
            span_b, span_b_row = y
            if are_overlapping(span_a, span_b):
                df.drop(span_b_row, axis=1, inplace=True)

    return df


def apply_annotations(writer, df):
    '''
        This takes a BioCWriter for 1 document and a Pandas Dataframe of annotations
        to apply to the BioCWriter

        Enforces as little DF modification as possible. This function is only
        intended to take the DF it was given and return a BioC Writer
    '''

    # If nothing in DF, we can safely return the non-modified writer
    # if df.shape[0] == 0 or df.shape[1] != 11:
    if df.shape[0] == 0:
        return writer

    # (TODO) Cases in which a user hasn't annotated something in the title...
    section_ids = list(df['section_id'].unique())
    if not len(section_ids) >= 1:
        raise ValueError('Incorrect number of document sections.')

    bioc_doc = writer.collection.documents[0]

    i = 0
    for offset_position, group_df in df.groupby(['section_offset']):
        # This is another approach to offset grouping:
        # offset = int(bioc_passage.offset)
        # for idx, row in df[df['offset'] == offset].iterrows():

        # (TODO) (WARNING) If only 1 section, then the bioc_passage will be wrong
        bioc_passage = bioc_doc.passages[i]
        i = i + 1

        # Should already be empty if from doc.as_writer(), but perhaps we add an
        # append method in the future
        bioc_passage.clear_annotations()

        for row_idx, row in group_df.iterrows():
            annotation = BioCAnnotation()

            annotation.id = str(row_idx)

            annotation.put_infon('uid', row['uid'])
            annotation.put_infon('source', row['source'])
            annotation.put_infon('user_id', str(row['user_id']))
            annotation.put_infon('type', row['ann_type'])
            annotation.put_infon('type_id', str(row['ann_type_id']))

            location = BioCLocation()
            location.offset = str(int(row['start_position']))
            location.length = str(int(row['length']))
            annotation.add_location(location)

            annotation.text = row['text']

            bioc_passage.add_annotation(annotation)

    return writer

