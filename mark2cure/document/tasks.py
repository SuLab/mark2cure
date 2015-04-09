from django.conf import settings

from mark2cure.document.models import Document, Section
from mark2cure.common.formatter import pad_split

from Bio import Entrez, Medline
from celery import task

import requests
from datetime import datetime, timedelta
import time

@task()
def get_pubtator_response(doc_pk, api_ann, session_id, data, payload, idx):
    attempts = 5

    while idx < attempts:
        idx += 1
        time.sleep(idx * 10)
        results = requests.post(
            'http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/{session_id}/Receive/'.format(session_id=session_id),
            data=data, params=payload)
        print api_ann +' : '+ results.content

        if results.content != 'Not yet':
            idx = 100

    if results.content and idx == 100:
        doc = Document.objects.get(pk=doc_pk)
        if api_ann == 'tmChem':
            doc.pubtator_chem = results.content
        if api_ann == 'DNorm':
            doc.pubtator_disease = results.content
        if api_ann == 'GNormPlus':
            doc.pubtator_gene = results.content
        doc.save()


@task()
def get_pubmed_document(pubmed_id):
    Entrez.email = settings.ENTREZ_EMAIL
    h = Entrez.efetch(db='pubmed', id=[str(pubmed_id),], rettype='medline', retmode='text')
    records = Medline.parse(h)

    # Reference to abbreviations: http://www.nlm.nih.gov/bsd/mms/medlineelements.html
    for record in records:
        if record.get('TI') and record.get('AB') and record.get('PMID') and record.get('CRDT'):
            if Document.objects.pubmed_count(record.get('PMID')) is 0:
                title = ' '.join( pad_split(record.get('TI')) )
                abstract = ' '.join( pad_split(record.get('AB')) )

                doc, doc_c = Document.objects.get_or_create(document_id=record.get('PMID'))
                doc.title = title
                doc.source = 'pubmed'
                doc.save()

                sec, sec_c = Section.objects.get_or_create(kind='t', document=doc)
                sec.text = title
                sec.save()

                sec, sec_c = Section.objects.get_or_create(kind='a', document=doc)
                sec.text = abstract
                sec.save()

                doc.save(force_update=True)



@task
def get_pubmed_documents(terms=settings.ENTREZ_TERMS):
    Entrez.email = settings.ENTREZ_EMAIL

    for term in terms:
        h = Entrez.esearch(db='pubmed', retmax=settings.ENTREZ_MAX_COUNT, term=term)
        result = Entrez.read(h)
        ids = result['IdList']
        h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
        records = Medline.parse(h)

        # Reference to abbreviations: http://www.nlm.nih.gov/bsd/mms/medlineelements.html
        for record in records:
            if record.get('TI') and record.get('AB') and record.get('PMID') and record.get('CRDT'):
                if Document.objects.pubmed_count(record.get('PMID')) is 0:
                    doc = Document.objects.create(document_id=record.get('PMID'))
                    doc.title = record.get('TI')
                    doc.created = datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
                    doc.source = 'pubmed'
                    doc.save()

                    sec = Section(kind='o')
                    sec.document = doc
                    sec.save()

                    sec = Section(kind='t')
                    sec.text = record.get('TI')
                    sec.document = doc
                    sec.save()

                    sec = Section(kind='a')
                    sec.text = record.get('AB')
                    sec.document = doc
                    sec.save()

