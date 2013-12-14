from django.conf import settings
from datetime import datetime
from django.core.exceptions import ObjectDoesNotExist
from django.contrib.auth.models import User

from mark2cure.document.models import Document, Section, View, Annotation
from Bio import Entrez, Medline
from bs4 import BeautifulSoup, NavigableString

import csv, random, re, nltk, datetime

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


def import_golden_documents():
    path = "NCBI_corpus_development"
    with open('assets/NCBI_corpus/'+ path +'.txt','r') as f:
        reader = csv.reader(f, delimiter='\t')
        for num, title, text in reader:
            try:
                doc = Document.objects.get(document_id = num)
            except ObjectDoesNotExist:
                doc = Document()

                doc.document_id = num
                doc.title = title
                doc.source = path
                doc.save()

                sec = Section(kind = "t")
                sec.text = title
                sec.document = doc
                sec.save()

                sec = Section(kind = "a")
                sec.text = text
                sec.document = doc
                sec.save()


def annotate_golden_documents():
    user, created = User.objects.get_or_create(username="goldenmaster")
    if created:
        user.set_password('')
        user.save()

    path = "NCBI_corpus_development"
    with open('assets/NCBI_corpus/'+ path +'_annos.txt','rU') as f:
        reader = csv.reader(f, delimiter='\t')
        for doc_id, doc_field, kind, text, start, stop in reader:
            doc = Document.objects.get(document_id = doc_id)
            for section in doc.section_set.all():
                if section.kind == doc_field[0]:
                    view, created = View.objects.get_or_create(section = section, user = user)

                    ann, created = Annotation.objects.get_or_create(view = view, text = text, start = start, type = kind)
                    ann.kind = "e"
                    ann.user_agent = "goldenmaster"
                    ann.save()


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


def strip_tags(html, invalid_tags):
    soup = BeautifulSoup(html)
    # print soup

    for tag in soup.findAll(True):
        if tag.name in invalid_tags:
            s = ""

            for c in tag.contents:
                if not isinstance(c, NavigableString):
                    c = strip_tags(unicode(c), invalid_tags)
                s += unicode(c)

            tag.replaceWith(s)

    return soup

