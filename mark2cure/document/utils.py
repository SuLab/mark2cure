from django.conf import settings
from datetime import datetime

from mark2cure.document.models import Document
from Bio import Entrez, Medline
from bs4 import BeautifulSoup, NavigableString

import random, re, nltk, datetime

def import_document():
    Document.objects.create_from_pubmed_id('24260557')


def get_pubmed_document(pubmed_id):
    Entrez.email = settings.ENTREZ_EMAIL
    h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
    records = Medline.parse(h)
    # for record in records:
      # if record.get('TI') and record.get('AB') and record.get('PMID') and record.get('CRDT'):
    return records


def get_pubmed_documents(terms = settings.ENTREZ_TERMS):
    Entrez.email = settings.ENTREZ_EMAIL

    for term in terms:
      h = Entrez.esearch(db='pubmed', retmax=settings.ENTREZ_MAX_COUNT, term=term)
      result = Entrez.read(h)
      ids = result['IdList']
      h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
      records = Medline.parse(h)

      #
      # Reference to abbreviations: http://www.nlm.nih.gov/bsd/mms/medlineelements.html
      #
      for record in records:
        if record.get('TI') and record.get('AB') and record.get('PMID') and record.get('CRDT'):
          print Document.objects.pubmed_count( record.get('PMID') )
          if Document.objects.pubmed_count( record.get('PMID') ) is 0:
            doc = Document.objects.create(document_id = record.get('PMID'))
            doc.title = record.get('TI')
            doc.text = record.get('AB')
            doc.created = datetime.datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
            doc.save()
