from django.conf import settings

from mark2cure.document.models import Document, Section

from Bio import Entrez, Medline
from celery import task
import datetime


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
                    doc.created = datetime.datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
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

