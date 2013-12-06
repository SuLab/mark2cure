from django.conf import settings
from datetime import datetime
from django.core.exceptions import ObjectDoesNotExist

from mark2cure.document.models import Document, Section
from Bio import Entrez, Medline
from bs4 import BeautifulSoup, NavigableString

import random, re, nltk, datetime

def create_from_pubmed_id(pubmed_id=None):
    pubmed_id = str(pubmed_id)

    ## Check if the account already exists
    try:
        doc = Document.objects.get(document_id = pubmed_id)
    except ObjectDoesNotExist:
        doc = Document()

        Entrez.email = settings.ENTREZ_EMAIL
        h = Entrez.efetch(db='pubmed', id=[pubmed_id], rettype='medline', retmode='text')
        records = Medline.parse(h)

        for record in records:
          if record.get('TI') and record.get('AB') and record.get('PMID') and record.get('CRDT'):
            doc.document_id = record.get('PMID')
            doc.title = record.get('TI')
            doc.created = datetime.datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
            doc.source = "pubmed"
            doc.save()

            sec = Section(kind = "t")
            sec.text = record.get('TI')
            sec.document = doc
            sec.save()

            sec = Section(kind = "a")
            sec.text = record.get('AB')
            sec.document = doc
            sec.save()
          break

    return doc


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
          if Document.objects.pubmed_count( record.get('PMID') ) is 0:
            doc = Document.objects.create(document_id = record.get('PMID'))
            doc.title = record.get('TI')
            doc.created = datetime.datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
            doc.source = "pubmed"
            doc.save()

            sec = Section(kind = "t")
            sec.text = record.get('TI')
            sec.document = doc
            sec.save()

            sec = Section(kind = "a")
            sec.text = record.get('AB')
            sec.document = doc
            sec.save()
