from mark2cure.common.bioc import BioCWriter, BioCCollection, BioCAnnotation, BioCLocation
from mark2cure.document.models import Document, Annotation

import datetime
import xmltodict
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
    text = text.replace("\\;", " ; ");
    text = text.replace("\\\"", " \" ")
    text = text.replace("  ", " ")
    text = text.replace("  ", " ")
    return nltk.word_tokenize( text.encode('utf-8') )


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
    return json.dumps(o)


def apply_bioc_annotations(writer, user=None):
    bioc_doc = writer.collection.documents[0]

    for bioc_passage in bioc_doc.passages:
        section_pk = bioc_passage.infons['id']
        offset = int(bioc_passage.offset)

        query_set = Annotation.objects.filter(view__section__pk=section_pk)
        if user:
            query_set = query_set.filter(view__user=user)

        for ann in query_set.all():
            annotation = BioCAnnotation()
            annotation.id = str(ann.pk)
            annotation.put_infon('user', str(ann.view.user.pk))
            annotation.put_infon('user_name', str(ann.view.user.username))

            # (TODO) Map type strings back to 0,1,2
            annotation.put_infon('type', str(ann.type))
            annotation.put_infon('type', str(0))

            location = BioCLocation()
            location.offset = str(offset + ann.start)
            location.length = str(len(ann.text))
            annotation.add_location(location)

            annotation.text = ann.text

            bioc_passage.add_annotation(annotation)

    return writer

