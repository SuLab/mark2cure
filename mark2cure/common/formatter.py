from django.contrib.contenttypes.models import ContentType

from ..common.bioc import BioCWriter, BioCCollection, BioCAnnotation, BioCLocation
from ..task.entity_recognition.models import EntityRecognitionAnnotation

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


def pubtator_df_as_writer(writer, df):
    bioc_doc = writer.collection.documents[0]
    approved_types = ['Disease', 'Gene', 'Chemical']

    i = 0
    for offset, group_df in df.groupby(['offset']):
        bioc_passage = bioc_doc.passages[i]
        i = i + 1

        bioc_passage.annotations = []

        group_df = clean_passage_df_overlaps(group_df)

        for row_idx, row in group_df.iterrows():

            if row['ann_type'] in approved_types:
                annotation = BioCAnnotation()
                annotation.id = str(row_idx)
                annotation.put_infon('user', 'pubtator')
                annotation.put_infon('uid', str(row['uid']))

                # (TODO) Map type strings back to 0,1,2
                annotation.put_infon('type', str(approved_types.index(row['ann_type'])))

                location = BioCLocation()
                loc = row['location'].split(':')
                location.offset = str(loc[0])
                location.length = str(loc[1])
                annotation.add_location(location)

                annotation.text = row['text']

                bioc_passage.add_annotation(annotation)

    return writer


def apply_bioc_annotations(writer, user=None):
    bioc_doc = writer.collection.documents[0]
    approved_types = ['disease', 'gene_protein', 'drug']

    for bioc_passage in bioc_doc.passages:
        section_pk = bioc_passage.infons['id']
        offset = int(bioc_passage.offset)
        content_type_id = str(ContentType.objects.get_for_model(EntityRecognitionAnnotation.objects.all().first()).id)
        if user:
            er_ann_query_set = EntityRecognitionAnnotation.objects.annotations_for_section_and_user(section_pk, user.pk, content_type_id)
        else:
            er_ann_query_set = EntityRecognitionAnnotation.objects.annotations_for_section_pk(section_pk, content_type_id)

        for er_ann in er_ann_query_set:
            annotation = BioCAnnotation()
            annotation.id = str(er_ann.id)
            annotation.put_infon('user', str(er_ann.user_id))

            # (TODO) Map type strings back to 0,1,2
            annotation.put_infon('type', str(approved_types.index(er_ann.type)))

            location = BioCLocation()
            location.offset = str(offset + er_ann.start)
            location.length = str(len(er_ann.text))
            annotation.add_location(location)

            annotation.text = er_ann.text

            bioc_passage.add_annotation(annotation)

    return writer

