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

        if len(user_annotations) is 0:
            raise ValueError("No user annotations available for this document")

        if len(gold_annotations) is 0:
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

      # print matches
      # print misses
      # print returned
      # print all_returned

      precision = matches / float(returned)
      recall = matches / float(all_returned)

      # print int(precision)
      # print recall

      if int(precision+recall) is not 0:
        f = ( 2 * precision * recall ) / ( precision + recall )
        return f
      else:
        pass

    def all_for_user(self, user_id):
      documents = db.session.query(Document).filter_by(source = 'NCBI_corpus_development').all()
      user = db.session.query(User).get( int(user_id) )
      results = []

      for doc in documents:
        user_annotations = db.session.query(Annotation).filter_by(document = doc).filter_by(user = user).all()
        if len(user_annotations) is not 0:
          results.append( self.user_vs_gold(user, doc) )


      print results
      results = [x for x in results if x is not None]
      print sum(results) / float(len(results))


    def run(self, document, user):
      self.all_for_user(user)

      # document = db.session.query(Document).get( int(document) )
      # user = db.session.query(User).get( int(user) )
      # print self.user_vs_gold(user, document)
