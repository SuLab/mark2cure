# -*- coding: utf-8 -*-
"""
    mark2cure.manage.analysis
    ~~~~~~~~~~~~~~~~~~~~~

    results/analysis commands
"""

from flask.ext.script import Command, Option, prompt, prompt_pass
from mark2cure.settings import *

from ..core import db
from ..models import User, Document, Annotation, View

import requests, re

def gold_matches(current_user, document):
    user_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user = current_user).all()
    gold_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user_id = 2).all()

    user_annotations = [ann.compare_view() for ann in user_annotations]
    gold_annotations = [ann.compare_view() for ann in gold_annotations]

    matches = len([ann for ann in user_annotations if ann in gold_annotations])

    return len(matches)

class Compare(Command):
    "F Score"
    option_list = (
        Option('--document', '-d', dest='document'),
        Option('--user', '-u', dest='user'),
    )

    def user_vs_gold(self, user, doc):
        # This is a document that requires validation
        user_annotations = db.session.query(Annotation).filter_by(document = doc).filter_by(user = user).all()
        gold_annotations = db.session.query(Annotation).filter_by(document = doc).filter_by(user_id = 2).all()

        user_annotations = [ann.compare_view() for ann in user_annotations]
        gold_annotations = [ann.compare_view() for ann in gold_annotations]

        if len(gold_annotations):
            print user_annotations
            print "/ / / / / / / / / "
            print gold_annotations
        else:
            raise ValueError("No golden annotations available for this document")

        return self.calc_f(user_annotations, gold_annotations)

    def calc_f(self, annotations_a, annotations_b):
      # It considers both the precision p and the recall r of the test to compute the score:
      # p is the number of correct results divided by the number of all returned results
      # r is the number of correct results divided by the number of results that should have been returned.
      # The F1 score can be interpreted as a weighted average of the precision and recall, where an F1 score reaches its best value at 1 and worst score at 0.
      matches = len([ann for ann in annotations_a if ann in annotations_b])
      misses = len(annotations_a) - matches

      # The count of the # of User's annotations
      returned = len(annotations_a)
      # The count of the # of GM annotations
      all_returned = len(annotations_b)

      # (p, r)
      return (matches/returned, matches/all_returned)

    def all(self):
      documents = db.session.query(Document).filter_by(source = 'NCBI_corpus_development').all()
      for doc in documents:
        gold_annotations = db.session.query(Annotation).filter_by(document = doc).filter_by(user_id = 1).all()
        print len(gold_annotations)



    def run(self, document, user):
      document = db.session.query(Document).get( int(document) )
      user = db.session.query(User).get( int(user) )

      print self.user_vs_gold(user, document)
