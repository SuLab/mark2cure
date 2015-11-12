from django.shortcuts import get_object_or_404
from django.http import HttpResponse

from mark2cure.common.formatter import bioc_writer, bioc_as_json
from .models import Document, Pubtator


def read_bioc(request, pubmed_id, format_type):
    """A plan (no annotations) BioC file for a PMID"""
    writer = bioc_writer(request)
    doc = get_object_or_404(Document, document_id=pubmed_id)

    writer = bioc_writer(request)
    bioc_document = doc.as_bioc_with_passages()
    writer.collection.add_document(bioc_document)

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')


def read_pubtator(request, pk):
    """Return the exact Pubtator file that was downloaded"""
    pubtator = get_object_or_404(Pubtator, pk=pk)
    return HttpResponse(pubtator.content, content_type='text/xml')


def read_pubtator_bioc(request, pubmed_id, format_type):
    """A merged file of the multiple Pubtator responses"""
    # When fetching via pubmed, include no annotaitons
    doc = get_object_or_404(Document, document_id=pubmed_id)

    writer = bioc_writer(request)
    bioc_document = doc.as_bioc_with_pubtator_annotations()
    writer.collection.add_document(bioc_document)

    if format_type == 'json':
        writer_json = bioc_as_json(writer)
        return HttpResponse(writer_json, content_type='application/json')
    else:
        return HttpResponse(writer, content_type='text/xml')

