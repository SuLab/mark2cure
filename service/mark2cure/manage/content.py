# -*- coding: utf-8 -*-
"""
    overholt.manage.users
    ~~~~~~~~~~~~~~~~~~~~~

    user management commands
"""
from flask.ext.script import Command, Option, prompt, prompt_pass
from mark2cure.settings import *
from ..models import *
from ..core import db

import requests, re, csv, datetime
import xml.etree.ElementTree as ET
# from Bio import Entrez, Medline
from bs4 import BeautifulSoup, NavigableString

# d = feedparser.parse('http://www.ncbi.nlm.nih.gov/entrez/eutils/erss.cgi?rss_guid=1ZIrvZUsS72fd3vWhvxLcpvS7DeyhJ23q0ooCoFlM1kM4TxkAL')
# pubmed_ids = []
# for paper in d.entries:
#   pubmed_ids.append( re.sub(r'\D', '', paper['id']) )
#
# Entrez.email = 'nanis@scripps.edu'
# records = Medline.parse( Entrez.efetch(db='pubmed', id=pubmed_ids, rettype='medline', retmode='text') )
#
# for record in records:
#     if record.get('TI') and record.get('AB') and record.get('PMID'):
#         print record.get('TI') # title
#         print record.get('AB') # abstract
#         print record.get('PMID') # pubmed id

class Import(Command):

    def run(self):

        pass
#     MAX_COUNT = 100
#     Entrez.email = 'nanis@scripps.edu'
#
#     for term in ['chordoma', 'cancer']:
#       h = Entrez.esearch(db='pubmed', retmax=MAX_COUNT, term=term)
#       result = Entrez.read(h)
#       ids = result['IdList']
#       h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
#       records = Medline.parse(h)
#
#       #
#       # Reference to abbreviations: http://www.nlm.nih.gov/bsd/mms/medlineelements.html
#       #
#       for record in records:
#         if record.get('TI') and record.get('AB') and record.get('PMID') and record.get('CRDT'):
#           instance = db.session.query(Document).filter_by(document_id = record.get('PMID')).count()
#           if instance == 0:
#             doc = Document( record.get('PMID'),
#                             record.get('AB'),
#                             record.get('TI'),
#                             datetime.datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
#                           )
#             db.session.add(doc)
#           db.session.commit()
    # print "Import"

def strip_tags(html, invalid_tags):
    soup = BeautifulSoup(html)

    for tag in soup.findAll(True):
        if tag.name in invalid_tags:
            s = ""

            for c in tag.contents:
                if not isinstance(c, NavigableString):
                    c = strip_tags(unicode(c), invalid_tags)
                s += unicode(c)

            tag.replaceWith(s)

    return soup

class SolidGold(Command):
    "Populate the documents with the NCBI SolidGold"

    option_list = (
        Option('--path', '-p', dest='path'),
    )

    def run(self, path):
      # Get the gold bot account to make the db entries
      user = db.session.query(User).get(2)
      with open('../assets/NCBI_corpus/'+ path +'.txt','r') as f:
        reader = csv.reader(f, delimiter='\t')
        for num, title, text in reader:
          text = text.replace('<category=', '<category type=')
          plain_text = str(strip_tags(text, ['category']))

          title = title.replace('<category=', '<category type=')
          plain_title = str(strip_tags(title, ['category']))

          soup = BeautifulSoup(text)

          # Add the document to the database
          doc = Document( int(num),
                          plain_text,
                          plain_title,
                          datetime.datetime.utcnow(),
                          path )
          db.session.add(doc)
          db.session.commit()

          cats = set([cat.text for cat in soup.findAll('category')])
          for cat in cats:
            for m in re.finditer( cat, plain_text ):
              # print( cat.text, cat['type'], m.start(), m.end() )
              ann = Annotation( 0,
                                'disease',
                                cat,
                                m.start(),
                                len(cat),
                                m.end(),
                                user,
                                doc,
                                'gold_1.0',
                                '',
                                None );
              db.session.add(ann)
              # Save every document instead of once incase some doc crashes
              db.session.commit()
