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
    '''
      Used on the document API to check for good performance from users
    '''
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

    def process_annotations(self, user=None, document=None, view=None):
      if user is not None and document is not None:
        annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user = user).all()
      elif view is not None:
        annotations = db.session.query(Annotation).filter_by(document = view.document).filter_by(user = view.user).all()
      else:
        raise ValueError("Not enough information given to retrieve annotations")

      # if len(annotations) is 0:
        # raise ValueError( "No user ({}) annotations available for this document ({})".format(user.username, document.document_id) )

      annotations = [ann.compare_view() for ann in annotations]
      # (TODO) FIGURE THIS OUT
      # Return back a uniq list of dictionaries
      return [dict(y) for y in set(tuple(x.items()) for x in annotations)]

    def exact(self, gm_ann, user_anns):
        '''
          Exact or Coextensive match finding for annotations.

          Works off start of annotation and cleaned length both being equal

        '''

        gm_end = gm_ann['start'] + len( gm_ann['text'].strip() )
        for user_ann in user_anns:
          if gm_ann['start'] is user_ann['start'] and gm_end is user_ann['start'] + len( user_ann['text'].strip() ):
            return True
        return False

    def partial(self, needle_ann, haystack_anns):
        # print needle_ann
        # print haystack_anns
        # print ""

        for compare_ann in haystack_anns:
          if distance.edit_distance( needle_ann.get('text'), compare_ann.get('text') ) <= 2:
            # If the user annotation is in one of the
            return True
          elif needle_ann.get('text') in compare_ann.get('text'):
            return True
        return False

    def calc_score(self, annotations_a, annotations_b, algo=0):
      '''
        This calculates the comparsion overlap between two arrays of dictionary terms

        Algorithms
          - 0: Exact
          - 1: Partial

        It considers both the precision p and the recall r of the test to compute the score:
        p is the number of correct results divided by the number of all returned results
        r is the number of correct results divided by the number of results that should have been returned.
        The F1 score can be interpreted as a weighted average of the precision and recall, where an F1 score reaches its best value at 1 and worst score at 0.

       tp  fp
       fn  tn

      '''

      # print "Theirs: {}".format( [( x.get('text'), x.get('start') ) for x in annotations_a] )
      # print "GM: {}".format( [( x.get('text'), x.get('start') ) for x in annotations_b] )

      # Correct annotations the user submitted
      if algo is 0:
        # For each of the GMs, check to see if the user has a coextensive match
        true_positive = [gm_ann for gm_ann in annotations_b if self.exact(gm_ann, annotations_a)]
      elif algo is 1:
        true_positive = [ann for ann in annotations_a if self.partial(ann, annotations_b)]
      else:
        true_positive = [gm_ann for gm_ann in annotations_b if self.exact(gm_ann, annotations_a)]

      true_positive = min(len(annotations_b), len(true_positive))

      # Annotations the user submitted that were wrong
      false_positive = len(annotations_a) - true_positive

      # Because of the potential for more true positive than false neg we need to
      # limit to 0, could also cap true-pos count of gm annotation
      false_negative = max(0, len(annotations_b) - true_positive )

      # print "TP: {} FP: {} FN: {}  || User Ann: {} GM Ann: {}".format(true_positive, false_positive, false_negative, len(annotations_a), len(annotations_b))
      return (true_positive, false_positive, false_negative)

    def determine_f(self, true_positive, false_positive, false_negative):
      if true_positive + false_positive is 0:
        return (0,0,0)

      precision = true_positive / float(true_positive + false_positive)
      recall = true_positive / float(true_positive + false_negative)

      # print "Precision: {} Recall: {}".format(precision, recall)

      if int(precision+recall) is not 0:
        f = ( 2 * precision * recall ) / ( precision + recall )
        return (precision, recall, f)
      else:
        return (0,0,0)

    def figure_one(self, document):
      print "Document {}\t Exact \t\t Partial\n".format(document.document_id)
      print "\t".join(['user          ', 'p', 'r', 'f', 'p', 'r', 'f'])

      gm_annotations = self.process_annotations(user=User.query.get(2), document=document)
      ncbo_annotations = self.process_annotations(user=User.query.get(1), document=document)

      exact_truth_table = self.calc_score(ncbo_annotations, gm_annotations, 0)
      exact_scores = self.determine_f( exact_truth_table[0], exact_truth_table[1], exact_truth_table[2] )

      partial_truth_table = self.calc_score(ncbo_annotations, gm_annotations, 1)
      partial_scores = self.determine_f( partial_truth_table[0], partial_truth_table[1], partial_truth_table[2] )
      print "\t".join(["NCBO          ", "%.2f"%exact_scores[0], "%.2f"%exact_scores[1], "%.2f"%exact_scores[2], "%.2f"%partial_scores[0], "%.2f"%partial_scores[1], "%.2f"%partial_scores[2]])


      worker_views = db.session.query(View).filter( View.created >= '2013-10-27' ).filter( View.user.has(mturk=1) ).filter_by( document = document ).all()
      for k in worker_views:
        annotations = self.process_annotations(view=k)
        exact_truth_table = self.calc_score(annotations, gm_annotations, 0)
        exact_scores = self.determine_f( exact_truth_table[0], exact_truth_table[1], exact_truth_table[2] )

        partial_truth_table = self.calc_score(annotations, gm_annotations, 1)
        partial_scores = self.determine_f( partial_truth_table[0], partial_truth_table[1], partial_truth_table[2] )
        print "\t".join([k.user.username, "%.2f"%exact_scores[0], "%.2f"%exact_scores[1], "%.2f"%exact_scores[2], "%.2f"%partial_scores[0], "%.2f"%partial_scores[1], "%.2f"%partial_scores[2]])

      print "--------------------------------------------------------\n"


    def figure_two(self):
      documents = db.session.query(Document).filter_by(source = 'NCBI_corpus_development').all()

      ncbo_list = []
      ncbo_part_list = []
      mturk_list = []
      mturk_part_list = []

      for document in documents:
        gm_annotations = self.process_annotations(user=User.query.get(2), document=document)
        ncbo_annotations = self.process_annotations(user=User.query.get(1), document=document)

        ncbo_list.append(       self.calc_score(ncbo_annotations, gm_annotations, 0) )
        ncbo_part_list.append(  self.calc_score(ncbo_annotations, gm_annotations, 1) )

        worker_views = db.session.query(View).filter( View.created >= '2013-10-27' ).filter( View.user.has(mturk=1) ).filter_by( document = document ).all()
        for k in worker_views:
          annotations = self.process_annotations(view=k)
          mturk_list.append(      self.calc_score(annotations, gm_annotations, 0) )
          mturk_part_list.append( self.calc_score(annotations, gm_annotations, 1) )


      ncbo_list = map(sum,zip(*ncbo_list))
      ncbo_part_list = map(sum,zip(*ncbo_part_list))
      mturk_list = map(sum,zip(*mturk_list))
      mturk_part_list = map(sum,zip(*mturk_part_list))

      ncbo_list = self.determine_f( ncbo_list[0], ncbo_list[1], ncbo_list[2] )
      ncbo_part_list = self.determine_f( ncbo_part_list[0], ncbo_part_list[1], ncbo_part_list[2] )
      mturk_list = self.determine_f( mturk_list[0], mturk_list[1], mturk_list[2] )
      mturk_part_list = self.determine_f( mturk_part_list[0], mturk_part_list[1], mturk_part_list[2] )

      # ncbo_list = [x/len(ncbo_list) for x in map(sum,zip(*ncbo_list))]
      # ncbo_part_list = [x/len(ncbo_part_list) for x in map(sum,zip(*ncbo_part_list))]
      # mturk_list = [x/len(mturk_list) for x in map(sum,zip(*mturk_list))]
      # mturk_part_list = [x/len(mturk_part_list) for x in map(sum,zip(*mturk_part_list))]
      # self.determine_f()

      print "--------------------------------------------------------\n"
      print "\t".join(['user ', 'p', 'r', 'f', 'p', 'r', 'f'])
      print "\t".join(["NCBO ", "%.2f"%ncbo_list[0], "%.2f"%ncbo_list[1], "%.2f"%ncbo_list[2], "%.2f"%ncbo_part_list[0], "%.2f"%ncbo_part_list[1], "%.2f"%ncbo_part_list[2]])
      print "\t".join(["MTURK", "%.2f"%mturk_list[0], "%.2f"%mturk_list[1], "%.2f"%mturk_list[2], "%.2f"%mturk_part_list[0], "%.2f"%mturk_part_list[1], "%.2f"%mturk_part_list[2]])
      print "--------------------------------------------------------\n"

    def select_worker_results(self):
      worker_views = db.session.query(View).filter( View.created >= '2013-10-27' ).filter( View.user.has(mturk=1) ).all()
      return len(worker_views)

    def run(self, document, user):
      # [self.figure_one(document) for document in db.session.query(Document).filter_by(source = 'NCBI_corpus_development').all()]
      # self.figure_two()
      print self.select_worker_results()
      return ":)"

    ############################################################
    # def all_for_user(self, user_id):
    #   documents = db.session.query(Document).filter_by(source = 'NCBI_corpus_development').all()
    #   user = db.session.query(User).get( int(user_id) )
    #   results = []

    #   for doc in documents:
    #     user_annotations = db.session.query(Annotation).filter_by(document = doc).filter_by(user = user).all()
    #     if len(user_annotations) is not 0:
    #       results.append( self.user_vs_gold(user, doc) )

    #   print results
    #   results = [x for x in results if x is not None]
    #   print sum(results) / float(len(results))

    # def all_mturk_gm_subs(self):
    #   gm_views = db.session.query(View).filter( View.created >= '2013-10-27' ).filter( View.user.has(mturk=1) ).filter( View.document.has(source='NCBI_corpus_development') ).all()
    #   gm_views = [view for view in gm_views if view.user.id != 244]
    #   print "Total quality MTurk GM Submissions {}".format( len(gm_views) )
    #   print ""

    #   agg_p = []
    #   agg_r = []
    #   agg_f = []
    #   for view in gm_views:
    #     scores = self.user_vs_gold(view.user, view.document)
    #     agg_p.append( scores[0] )
    #     agg_r.append( scores[1] )
    #     agg_f.append( scores[2] )
    #     print "Worker: {} Document: {} F-Score: {}".format(view.user.username, view.document.id, scores[2])
    #     print ""

    #   a_p = sum(agg_p) / float(len(agg_p))
    #   a_r = sum(agg_r) / float(len(agg_r))
    #   a_f = sum(agg_f) / float(len(agg_f))

    #   print "Aggregated Precision: {} Recall {} F-Score: {}".format( a_p, a_r, a_f )

    # def user_vs_gold(self, user, doc):
    #     # This is a document that requires validation
    #     user_annotations = db.session.query(Annotation).filter_by(document = doc).filter_by(user = user).all()
    #     gold_annotations = db.session.query(Annotation).filter_by(document = doc).filter_by(user_id = 2).all()

    #     user_annotations = [ann.compare_view() for ann in user_annotations]
    #     user_annotations = [dict(y) for y in set(tuple(x.items()) for x in user_annotations)]

    #     gold_annotations = [ann.compare_view() for ann in gold_annotations]

    #     if len(user_annotations) is 0:
    #         raise ValueError("No user annotations available for this document")

    #     if len(gold_annotations) is 0:
    #         raise ValueError("No golden annotations available for this document")

    #     return self.calc_f(user_annotations, gold_annotations)
