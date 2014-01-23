from django.conf import settings
from datetime import datetime
from django.core.exceptions import ObjectDoesNotExist
from django.contrib.auth.models import User

from mark2cure.document.models import Document, Section, View, Annotation
from mark2cure.common.utils import Turk

from Bio import Entrez, Medline
from bs4 import BeautifulSoup, NavigableString

import re, nltk, datetime

def gold_matches(current_user, document):
    '''
      Used on the document API to check for good performance from users

      This is very "loose" but given how good turkers have been and how careful we want to be it's perfect
    '''
    user_annotations = Annotation.objects.filter(view__section__document = document, view__user = current_user).all()
    gold_user, created = User.objects.get_or_create(username="goldenmaster")
    gold_annotations = Annotation.objects.filter(view__section__document = document, view__user = gold_user).all()

    user_annotations = [ann.text for ann in user_annotations]
    gold_annotations = [ann.text for ann in gold_annotations]
    true_positives = [gm_ann for gm_ann in gold_annotations if gm_ann in user_annotations]

    return len(true_positives)


def check_validation_status(user, document, view=None):
    views = View.objects.filter(user = user, section__document = document ).all()
    if len(views) > 2:
      if user.profile.mturk:
        t = Turk()
        t.mtc.block_worker(user.username, "Attempted to submit same document multiple times.")
      return "Cannot submit a document twice"
      raise ValueError("Cannot submit a document twice")

    if document.section_set.all()[:1].get().validate and user.profile.mturk:
      # If this is a validate document, check the user's history, if it's their 3rd submission
      # or more run test to potentially fail if poor performance
      if view is None:
        valid_views = View.objects.filter(user = user, section__validate = True).order_by('-created').all()[:6]
      else:
        valid_views = View.objects.filter(user = user, section__validate = True, created__lt = view.updated).order_by('-created').all()[:6]

      if len(valid_views) is 6:
        if sum(1 for x in valid_views if gold_matches(x.user, x.section.document) >= 1) < 3:
          t = Turk()
          t.mtc.block_worker(user.username, "Failed to properly answer golden master performance documents")
          return "Failed"
        else:
          return "Passed w/ ", sum(1 for x in valid_views if gold_matches(x.user, x.section.document) >= 1)
      else:
        return "Not enough yet"


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

