from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist
from django.contrib.auth.models import User
from django.db.models import Q

from mark2cure.document.models import Document, Section, View, Annotation

from datetime import datetime
from Bio import Entrez, Medline

import re, nltk, datetime


def determine_f(true_positive, false_positive, false_negative):
    if float(true_positive + false_positive) == 0.0:
      return (0.0, 0.0, 0.0)

    if float(true_positive + false_negative) == 0.0:
      return (0.0, 0.0, 0.0)

    precision = true_positive / float(true_positive + false_positive)
    recall = true_positive / float(true_positive + false_negative)

    if float(precision + recall) > 0.0:
      f = ( 2 * precision * recall ) / ( precision + recall )
      return (precision, recall, f)
    else:
      return (0.0, 0.0, 0.0)


def match_exact(gm_ann, user_anns):
    for user_ann in user_anns:
        if user_ann.is_exact_match(gm_ann): return True
    return False


def generate_results(document, user):
    '''
      This calculates the comparsion overlap between two arrays of dictionary terms

      It considers both the precision p and the recall r of the test to compute the score:
      p is the number of correct results divided by the number of all returned results
      r is the number of correct results divided by the number of results that should have been returned.
      The F1 score can be interpreted as a weighted average of the precision and recall, where an F1 score reaches its best value at 1 and worst score at 0.

     tp  fp
     fn  *tn

    '''
    gm_annotations = document.latest_annotations()
    user_annotations = document.latest_annotations(user)

    true_positives = [gm_ann for gm_ann in gm_annotations if match_exact(gm_ann, user_annotations)]

    # Annotations the user submitted that were wrong (the User set without their True Positives)
    # false_positives = user_annotations - true_positives
    false_positives = user_annotations
    for tp in true_positives:
      false_positives = false_positives.exclude(start = tp.start, text = tp.text)


    # # Annotations the user missed (the GM set without their True Positives)
    # false_negatives = gm_annotations - true_positives
    false_negatives = gm_annotations
    for tp in true_positives:
      false_negatives = false_negatives.exclude(start = tp.start, text = tp.text)

    score = determine_f( len(true_positives), false_positives.count(), false_negatives.count() )
    return ( score, true_positives, false_positives, false_negatives )


def create_from_pubmed_id(pubmed_id=None):
    pubmed_id = str(pubmed_id)

    ## Check if the account already exists
    try:
        return Document.objects.get(document_id = pubmed_id)
    except ObjectDoesNotExist:
        Entrez.email = settings.ENTREZ_EMAIL
        h = Entrez.efetch(db='pubmed', id=[pubmed_id], rettype='medline', retmode='text')
        records = Medline.parse(h)

        for record in records:
          # http://www.nlm.nih.gov/bsd/mms/medlineelements.html

          if record.get('TI') and record.get('PMID') and record.get('CRDT'):
            doc = Document()

            doc.document_id = record.get('PMID')
            doc.title = record.get('TI')
            doc.created = datetime.datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
            doc.source = "pubmed"
            doc.save()

            sec = Section(kind = "o")
            sec.document = doc
            sec.save()

            sec = Section(kind = "t")
            sec.text = record.get('TI')
            sec.document = doc
            sec.save()

            if record.get('AB'):
              sec = Section(kind = "a")
              sec.text = record.get('AB')
              sec.document = doc
              sec.save()

            return doc
          break


