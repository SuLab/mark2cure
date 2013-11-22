from django.conf import settings
from datetime import datetime

from Bio import Entrez, Medline
from bs4 import BeautifulSoup, NavigableString

import random, re, nltk

def get_pubmed_documents(terms = ['chordoma', 'cancer']):
    MAX_COUNT = 100
    Entrez.email = 'nanis@scripps.edu'

    for term in terms:
      h = Entrez.esearch(db='pubmed', retmax=MAX_COUNT, term=term)
      result = Entrez.read(h)
      ids = result['IdList']
      h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
      records = Medline.parse(h)

      #
      # Reference to abbreviations: http://www.nlm.nih.gov/bsd/mms/medlineelements.html
      #
      for record in records:
        if record.get('TI') and record.get('AB') and record.get('PMID') and record.get('CRDT'):
          # instance = db.session.query(Document).filter_by(document_id = record.get('PMID')).count()
          # if instance == 0:
            print record.get('PMID'), record.get('TI')
            print record.get('AB')
            print "\n\n - - - - - - - - - - - - \n\n"
            # doc = Document( record.get('PMID'),
            #                 record.get('AB'),
            #                 record.get('TI'),
            #                 datetime.datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
            #               )
