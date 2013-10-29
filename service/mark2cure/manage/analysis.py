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

from nltk.metrics import *

import requests, re

def gold_matches(current_user, document):
    user_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user = current_user).all()
    gold_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user_id = 2).all()

    user_annotations = [ann.compare_view() for ann in user_annotations]
    gold_annotations = [ann.compare_view() for ann in gold_annotations]

    return len([ann for ann in user_annotations if ann in gold_annotations])

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
        user_annotations = [dict(y) for y in set(tuple(x.items()) for x in user_annotations)]

        gold_annotations = [ann.compare_view() for ann in gold_annotations]

        if len(user_annotations) is 0:
            raise ValueError("No user annotations available for this document")

        if len(gold_annotations) is 0:
            raise ValueError("No golden annotations available for this document")

        return self.calc_f(user_annotations, gold_annotations)

    def compare(self, needle_ann, haystack_anns):
        # print needle_ann
        # print haystack_anns
        # print ""

        for compare_ann in haystack_anns:
          if distance.edit_distance( needle_ann.get('text'), compare_ann.get('text') ) <= 2:
            return True
          elif needle_ann.get('text') in compare_ann.get('text'):
            return True
        return False

    def calc_f(self, annotations_a, annotations_b):
      # It considers both the precision p and the recall r of the test to compute the score:
      # p is the number of correct results divided by the number of all returned results
      # r is the number of correct results divided by the number of results that should have been returned.
      # The F1 score can be interpreted as a weighted average of the precision and recall, where an F1 score reaches its best value at 1 and worst score at 0.
      print "Theirs: {}".format( [( x.get('text'), x.get('start') ) for x in annotations_a] )
      print "GM: {}".format( [( x.get('text'), x.get('start') ) for x in annotations_b] )

      #
      #  tp  fp
      #  fn  tn
      #

      # Correct annotations the user submitted
      true_positive = len([ann for ann in annotations_a if self.compare(ann, annotations_b)])
      # true_positive = len([ann for ann in annotations_a if ann in annotations_b])

      # Annotations the user submitted that were wrong
      false_positive = len(annotations_a) - true_positive
      false_negative = len(annotations_b) - true_positive

      print "TP: {} FP: {} FN: {}  || User Ann: {} GM Ann: {}".format(true_positive, false_positive, false_negative, len(annotations_a), len(annotations_b))

      precision = true_positive / float(true_positive + false_positive)
      recall = true_positive / float(true_positive + false_negative)

      print "Precision: {} Recall: {}".format(precision, recall)

      if int(precision+recall) is not 0:
        f = ( 2 * precision * recall ) / ( precision + recall )
        return (precision, recall, f)
      else:
        return (0,0,0)

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

    def all_mturk_gm_subs(self):
      gm_views = db.session.query(View).filter( View.user.has(mturk=1) ).filter( View.document.has(source='NCBI_corpus_development') ).all()
      gm_views = [view for view in gm_views if view.user.id != 244]
      print "Total quality MTurk GM Submissions {}".format( len(gm_views) )
      print ""

      agg_p = []
      agg_r = []
      agg_f = []
      for view in gm_views:
        scores = self.user_vs_gold(view.user, view.document)
        agg_p.append( scores[0] )
        agg_r.append( scores[1] )
        agg_f.append( scores[2] )
        print "Worker: {} Document: {} F-Score: {}".format(view.user.username, view.document.id, scores[2])
        print ""

      a_p = sum(agg_p) / float(len(agg_p))
      a_r = sum(agg_r) / float(len(agg_r))
      a_f = sum(agg_f) / float(len(agg_f))

      print "Aggregated Precision: {} Recall {} F-Score: {}".format( a_p, a_r, a_f )

    def run(self, document, user):
      # self.all_for_user(user)
      self.all_mturk_gm_subs()

      # document = db.session.query(Document).get( int(document) )
      # user = db.session.query(User).get( int(user) )
      # print self.user_vs_gold(user, document)
