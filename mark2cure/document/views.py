from django.shortcuts import get_object_or_404
from django.http import HttpResponse

from ..common.formatter import bioc_as_json, clean_df, apply_annotations
from .models import Document, Pubtator

from rest_framework.decorators import api_view
from rest_framework.response import Response


@api_view(['GET'])
def read_document(request, document_pk):
    '''
        Return a JSON response with the generic document structure

        No annotations of any kind are included
    '''
    doc = get_object_or_404(Document, pk=document_pk)
    data = Document.objects.as_json(documents=[doc])
    return Response(data[0])


def read_pubtator_bioc(request, pubmed_id):
    """
        A merged file of the multiple Pubtator responses
        Prevents overlap of text as this is used by YPet
    """
    # When fetching via pubmed, include no annotaitons
    doc = get_object_or_404(Document, document_id=pubmed_id)

    df = Document.objects.ner_df(documents=[doc], include_pubtator=True)
    df = clean_df(df, overlap_protection=True)

    # writer = Document.objects.as_writer(documents=[doc])
    # writer = apply_annotations(writer, er_df=df)
    print(df)

    data = Document.objects.as_json(documents=[doc])
    return Response(data)


def read_pubtator(request, pk):
    """
        Return the exact Pubtator file that was downloaded
    """
    pubtator = get_object_or_404(Pubtator, pk=pk)
    return HttpResponse(pubtator.content, content_type='text/xml')


