from django.shortcuts import get_object_or_404
from django.http import HttpResponse

from ..common.formatter import bioc_as_json, clean_df, apply_annotations
from .models import Document, Pubtator


def read_bioc(request, pubmed_id, format_type):
    '''
        A plain (no annotations) BioC file for a PMID
    '''
    doc = get_object_or_404(Document, document_id=pubmed_id)

    writer = doc.as_writer()

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')


def read_pubtator(request, pk):
    '''
        Return the exact Pubtator file that was downloaded
    '''
    pubtator = get_object_or_404(Pubtator, pk=pk)
    return HttpResponse(pubtator.content, content_type='text/xml')


def read_pubtator_bioc(request, pubmed_id, format_type):
    '''
        A merged file of the multiple Pubtator responses
        Prevents overlap of text as this is used by YPet
    '''
    # When fetching via pubmed, include no annotaitons
    doc = get_object_or_404(Document, document_id=pubmed_id)

    df = doc.as_er_df_with_pubtator_annotations()
    df = clean_df(df, overlap_protection=True)

    writer = doc.as_writer()
    writer = apply_er_annotations(writer, df)

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')

