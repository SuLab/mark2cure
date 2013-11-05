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

import requests, re, collections

def gold_matches(current_user, document):
    '''
      Used on the document API to check for good performance from users
    '''
    user_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user = current_user).all()
    gold_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user_id = 2).all()

    user_annotations = [ann.compare_view() for ann in user_annotations]
    gold_annotations = [ann.compare_view() for ann in gold_annotations]

    return len([ann for ann in user_annotations if ann in gold_annotations])

def make_hashable(d):
    return (frozenset(x.iteritems()) for x in d)

def subtract(a, b):
    """ Remove the keys in b from a. """
    for k in b:
        if k in a:
            if isinstance(b[k], dict):
                subtract(a[k], b[k])
            else:
                del a[k]

class Compare(Command):
    "F Score"
    option_list = (
        Option('--document', '-d', dest='document'),
        Option('--user', '-u', dest='user'),
    )

    def __init__(self):
      self.error_aggreements = {'false_positives': [], 'false_negatives': []}

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
          Exact or Coextensive match finding for annotations. Works off start of annotation and cleaned length both being equal

          Returns True is any of the user annotations are equal to this GM Annotation

        '''
        gm_len = len(gm_ann['text'])
        for user_ann in user_anns:
          if gm_ann['start'] == user_ann['start']:
            if gm_len == len(user_ann['text']):
              return True

        return False

    def partial(self, gm_ann, user_anns):
        '''
          Overlap (Partial) match finding for annotations.

          Works off start of annotation and cleaned length enclosing a user annotation

        '''
        gm_len = len( gm_ann['text'])
        gm_start = gm_ann['start']

        for user_ann in user_anns:
          if gm_start <= user_ann['start'] <= gm_start + gm_len:
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
      # TP is the # of GM annotations that the user had as well
      if algo is 0:
        # For each of the GMs, check to see if the user has a coextensive match
        true_positives = [gm_ann for gm_ann in annotations_b if self.exact(gm_ann, annotations_a)]
      elif algo is 1:
        true_positives = [gm_ann for gm_ann in annotations_b if self.partial(gm_ann, annotations_a)]
      else:
        true_positives = [gm_ann for gm_ann in annotations_b if self.exact(gm_ann, annotations_a)]

      # In order to make our comparisons we need to do the conversion to a list of tuples
      true_positives = [tuple(sorted(item.items())) for item in true_positives]
      annotations_a = [tuple(sorted(item.items())) for item in annotations_a]
      annotations_b = [tuple(sorted(item.items())) for item in annotations_b]

      # Annotations the user submitted that were wrong
      false_positives = set(annotations_a) - set(true_positives)

      # Annotations the user missed
      false_negatives = set(annotations_b) - set(true_positives)

      for fp in false_positives:
        self.error_aggreements['false_positives'].append(fp[1][1])

      for fn in false_negatives:
        self.error_aggreements['false_negatives'].append(fn[1][1])

      # print "TP: {} \nFP: {} \nFN: {} \n|| \nUser Ann: {} \nGM Ann: {} \n\n-----\n\n".format(true_positive, false_positive, false_negative, annotations_a, annotations_b)
      return ( len(true_positives), len(false_positives), len(false_negatives) )

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

    def figure_two(self, documents, match):
      '''
        Generates a table of NCBO vs K1-5 users annotations for the list of input Documents

        Input:
          documents - List of SQLAlchemy Document models
          match - int for type of matching algorithmn to be used. 0 = Exact, 1 = Partial
        Output: none (prints ASCII Table)
      '''
      # Setup the results dictionary witch will store all the annotations we're comparing
      results = {'ncbo':[]}
      for item in range(0,6): results[item] = []

      for document in documents:
        # Collect the list of Annotations models for the Golden Master and NCBO Annotator to use throughout
        gm_annotations = self.process_annotations(user=User.query.get(2), document=document)
        ncbo_annotations = self.process_annotations(user=User.query.get(1), document=document)
        results['ncbo'].append( self.calc_score(ncbo_annotations, gm_annotations, match) )

        # Collect all of the MTurk workers that viewed this document
        worker_views = db.session.query(View).\
                        filter( View.created >= '2013-10-27' ).\
                        filter( View.user.has(mturk=1) ).\
                        filter_by( document = document ).\
                        all()

        workers_culmulative = []
        for worker_view in worker_views:
          annotations = self.process_annotations(view=worker_view)
          workers_culmulative.append( annotations )
        # Flatten the list
        workers_culmulative = [item for sublist in workers_culmulative for item in sublist]

        # looping through the uniq anns to get their counts
        b = {}
        for item in range(0,6):
          b[item] = []

        for ann in [dict(y) for y in set(tuple(x.items()) for x in workers_culmulative)]:
          for item in range(0, workers_culmulative.count(ann)):
            b[ item ].append( ann )

        # take the document results and put them into the global results
        for item in range(0,6):
          results[item].append( self.calc_score(b[item], gm_annotations, match) )

      print "--------------------------------------------------------\n"
      print "\t".join(['user ', 'p', 'r', 'f'])
      for i in results.iterkeys():
        results[i] = map(sum,zip(*results[i]))
        results[i] = self.determine_f( results[i][0], results[i][1], results[i][2] )
        print "\t".join(["{} ".format(i), "%.2f"%results[i][0], "%.2f"%results[i][1], "%.2f"%results[i][2]])
      print "--------------------------------------------------------\n"

    def select_worker_results(self):
      worker_views = db.session.query(View).\
          filter( View.created >= '2013-10-27' ).\
          filter( View.user.has(mturk=1) ).\
          all()
      return len(worker_views)

    def show_missed_results(self):
      '''
        Prints out the ordered list of the top False Positives and False Negatives
        that have been encountered globally when running the analysis
      '''
      print 'False_positives:'
      print collections.Counter(self.error_aggreements['false_positives'])
      print '\n\nFalse_negatives:'
      print collections.Counter(self.error_aggreements['false_negatives'])

    def run(self, document, user):
      documents = db.session.query(Document).\
          filter_by(source = 'NCBI_corpus_development').\
          limit(10).\
          all()
      self.figure_two(documents, 0)

      # self.show_missed_results()
      # print self.select_worker_results()
      return ":)"

    ############################################################

    # def figure_one(self, document):
    #   print "Document {}\t Exact \t\t Partial\n".format(document.document_id)
    #   print "\t".join(['user          ', 'p', 'r', 'f', 'p', 'r', 'f'])

    #   gm_annotations = self.process_annotations(user=User.query.get(2), document=document)
    #   ncbo_annotations = self.process_annotations(user=User.query.get(1), document=document)

    #   exact_truth_table = self.calc_score(ncbo_annotations, gm_annotations, 0)
    #   exact_scores = self.determine_f( exact_truth_table[0], exact_truth_table[1], exact_truth_table[2] )

    #   partial_truth_table = self.calc_score(ncbo_annotations, gm_annotations, 1)
    #   partial_scores = self.determine_f( partial_truth_table[0], partial_truth_table[1], partial_truth_table[2] )
    #   print "\t".join(["NCBO          ", "%.2f"%exact_scores[0], "%.2f"%exact_scores[1], "%.2f"%exact_scores[2], "%.2f"%partial_scores[0], "%.2f"%partial_scores[1], "%.2f"%partial_scores[2]])


    #   worker_views = db.session.query(View).filter( View.created >= '2013-10-27' ).filter( View.user.has(mturk=1) ).filter_by( document = document ).all()
    #   for k in worker_views:
    #     annotations = self.process_annotations(view=k)
    #     exact_truth_table = self.calc_score(annotations, gm_annotations, 0)
    #     exact_scores = self.determine_f( exact_truth_table[0], exact_truth_table[1], exact_truth_table[2] )

    #     partial_truth_table = self.calc_score(annotations, gm_annotations, 1)
    #     partial_scores = self.determine_f( partial_truth_table[0], partial_truth_table[1], partial_truth_table[2] )
    #     print "\t".join([k.user.username, "%.2f"%exact_scores[0], "%.2f"%exact_scores[1], "%.2f"%exact_scores[2], "%.2f"%partial_scores[0], "%.2f"%partial_scores[1], "%.2f"%partial_scores[2]])

    #   print "--------------------------------------------------------\n"

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
