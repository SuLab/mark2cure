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

class Compare(Command):
    "F Score"
    option_list = (
        Option('--document', '-d', dest='document'),
        Option('--user', '-u', dest='user'),
    )

    def __init__(self):
      self.error_aggreements = {'false_positives': [], 'false_negatives': []}

    def process_annotations(self, user=None, document=None, view=None, strict=False):
      '''
        This function returns all the dictionary annotations for a particular user and document or
        for a given View Model
      '''
      if user is not None and document is not None:
        annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user = user).all()
      elif view is not None:
        annotations = db.session.query(Annotation).filter_by(document = view.document).filter_by(user = view.user).all()
      else:
        raise ValueError("Not enough information given to retrieve annotations")

      if strict and len(annotations) is 0:
        raise ValueError( "No annotations available for this document")

      # Convert all of the Anotation SQLAlchemy models to cleaned up dictionaries with
      # cleaned whitespace and offsets
      annotations = [ann.compare_view() for ann in annotations]
      uniq_annotations = [dict(y) for y in set(tuple(x.items()) for x in annotations)]

      # if len(annotations) is not len(uniq_annotations):
      #   print annotations
      #   print uniq_annotations
      #   print ' - - - - - - - - - - '

      return uniq_annotations

    def exact(self, gm_ann, user_anns):
        '''
          Exact or Coextensive match finding for annotations. Works off start of annotation and cleaned length both being equal

          Returns True is any of the user annotations are equal to this GM Annotation

        '''
        gm_len = len(gm_ann['text'])
        for user_ann in user_anns:
          if gm_ann['start'] == user_ann['start'] and gm_len == len(user_ann['text']): return True
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
       fn  *tn

      '''
      # Calculate the correct annotations the user submitted TP is the # of GM
      # annotations that the user had as well
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

      # Annotations the user submitted that were wrong (the User set without their True Positives)
      false_positives = set(annotations_a) - set(true_positives)

      # Annotations the user missed (the GM set without their True Positives)
      false_negatives = set(annotations_b) - set(true_positives)

      # Add the False Positives and False Negatives to a global array to keep track which
      # annotations are commonly incorrect
      for fp in false_positives: self.error_aggreements['false_positives'].append(fp[1][1])
      for fn in false_negatives: self.error_aggreements['false_negatives'].append(fn[1][1])

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

        ncbo_score = self.calc_score(ncbo_annotations, gm_annotations, match)
        results['ncbo'].append( ncbo_score )

        # Collect all of the MTurk workers that viewed this document and assemble a list of all
        # the annotations across the MTurk Workers for this document in the experiment
        worker_views = db.session.query(View).\
                        filter( View.created >= '2013-10-27' ).\
                        filter( View.user.has(mturk=1) ).\
                        filter_by( document = document ).\
                        all()
        workers_culmulative = []
        for worker_view in worker_views:
          annotations = self.process_annotations(view=worker_view)
          workers_culmulative.append( annotations )
        # Flatten the list so we just have an array of all the combined workers annotations
        workers_culmulative = [item for sublist in workers_culmulative for item in sublist]

        # Build dictionary to store different K scores
        k = {}
        for item in range(0,6): k[item] = []
        # Looping through all the unique annotations to get their counts to actual
        # submitted annotations for the workers results
        for ann in [dict(y) for y in set(tuple(x.items()) for x in workers_culmulative)]:
          # Put that annotation into the k for the # that it matches (how many times did workers agree on that
          # particular annotation) and everything below it
          # Ex: If an annotation matches 3 times, it also matches 2 times
          for item in range(0, workers_culmulative.count(ann)): k[ item ].append( ann )

        # Put the document K score annotations and append their TP/FP/FN counts to the K results
        for i in range(0,6): results[i].append( self.calc_score(k[i], gm_annotations, match) )

      # We've now built up the results dictionary for our K scores and NCBO annotator for all the documents.
      # Sum all the scores up, calculate their P/R/F and print it out
      print "--------------------------------------------------------\n"
      print "\t".join(['user ', 'p', 'r', 'f'])
      for group in sorted(results.iterkeys()):
        results[group] = map(sum,zip(*results[group]))
        results[group] = self.determine_f( results[group][0], results[group][1], results[group][2] )
        print "\t".join(["{} ".format(group), "%.2f"%results[group][0], "%.2f"%results[group][1], "%.2f"%results[group][2]])
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
      # documents = db.session.query(Document).\
      #     filter_by(source = 'NCBI_corpus_development').\
      #     all()
      # self.figure_two(documents, 0)

      # self.show_missed_results()
      print self.select_worker_results()
      return ":)"
