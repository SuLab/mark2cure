from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
# from django.db.models import Count
from django.contrib.auth.models import User

from mark2cure.document.models import Document, Section, View, Annotation
from mark2cure.document.utils import check_validation_status
from mark2cure.account.models import Ncbo

from mark2cure.common.utils import Turk

import os, os.path, csv

class Command(BaseCommand):
    args = '<experiment_run_id>'
    help = 'Run and calculate metrics for the experiment'

    def handle(self, *args, **options):
        if len(args) < 2: raise Exception('Analysis needs 2 parameters <experiement_id, command>')

        experiment = args[0]
        command = args[1]

        self.stdout.write('-- Running analytics ({0}) on Experiment Run {1} --'.format(command, experiment))

        try:
          base_path = os.path.join(settings.PROJECT_PATH, 'results')
          os.chdir(base_path)
        except RuntimeError:
          print "A 'results' folder must be created in the project directory"

        if command == "ncbo_spectrum":
          documents = Document.objects.filter(source = 'NCBI_corpus_development').all()
          self.util_ncbo_specturm(documents)

        elif command == "worker_specturm":
          documents = Document.objects.filter(source = 'NCBI_corpus_development').all()
          self.util_worker_specturm(documents)

        elif command == "worker_ban_analysis":
          self.util_worker_ban_analysis()

        elif command == "worker_qualification_scores":
          self.util_worker_qualification_scores()

        else:
          pass


    def util_worker_qualification_scores(self):
        workers = User.objects.filter(userprofile__mturk = True).all()
        turk = Turk()
        results = []
        with open('mturk_qual_3_summary.tsv', 'wb') as csvfile:
          writer = csv.writer(csvfile, delimiter='\t')
          for page in range(1,5):
            quals = turk.mtc.get_qualifications_for_qualification_type(settings.AWS_QUAL_TEST_3, page_number = page)
            for qi in quals:
              writer.writerow([qi.Status, qi.SubjectId, qi.IntegerValue, qi.GrantTime])


    def util_worker_ban_analysis(self):
        workers = User.objects.filter(userprofile__mturk = True).all()
        for worker in workers:
          print "\n -- User "+ str(worker.id) +" -- \n"
          views = View.objects.filter(user = worker).all()
          for view in views:
            print view.section.validate, ' :: ', view.section.document.id, " :: ", check_validation_status(worker, view.section.document, view)



    def util_worker_specturm(self, documents):
        '''
          Generates a table of NCBO vs K1-5 users annotations for the list of input Documents
        '''
        gm_user = User.objects.get(username__exact = 'goldenmaster')
        types = ["*", "disease:modifier", "disease:class", "disease:specific", "disease:composite"]

        for disease_type in types:


          results = {}
          for i in range(1, 10):
            results[i] = {}
            results[i]['true_positives'] = []
            results[i]['false_positives'] = []
            results[i]['false_negatives'] = []


          for document in documents:
            for section in document.section_set.all():

              # Get the annotation for this grouping
              if disease_type is "*":
                gold_query = Annotation.objects.filter(view__section = section, view__user = gm_user)
                work_query = Annotation.objects.filter(view__section = section, experiment = 3, view__user__userprofile__mturk = True)
              else:
                gold_query = Annotation.objects.filter(type = disease_type, view__section = section, view__user = gm_user)
                work_query = Annotation.objects.filter(type = 'disease', view__section = section, experiment = 3, view__user__userprofile__mturk = True)


              # When gold_k is 0, that section (likely a title) has no disease terms in it
              gold_k = len(gold_query.values('view__user').distinct())
              worker_k = len(work_query.values('view__user').distinct())

              # Looping through all the unique annotations to get their counts to actual
              # submitted annotations for the workers results
              worker_dict_anns = list( work_query.values('text', 'start') )
              gm_dict_anns = list( gold_query.values('text', 'start') )

              worker_anns = work_query.all()
              worker_uniq_anns = worker_anns.values('text', 'start').distinct()

              k = {}
              for i in range(1, worker_k + 10): k[i] = []

              for ann in worker_uniq_anns:
                # Put that annotation into the k for the # that it matches (how many times did workers agree on that
                # particular annotation) and everything below it
                # Ex: If an annotation matches 3 times, it also matches 2 times
                for k_group in range(1, worker_dict_anns.count(ann) + 2 ): k[k_group].append( ann )


              # For this section, go back and actually store the results
              for k_group in range(1, worker_k + 1):
                # Put the document K score annotations and append their TP/FP/FN counts to the K results
                truths = self.calc_truth_table( k[k_group], gm_dict_anns)

                results[k_group]['true_positives'] = results[k_group]['true_positives'] + truths[0]
                results[k_group]['false_positives'] = results[k_group]['false_positives'] + truths[1]
                results[k_group]['false_negatives'] = results[k_group]['false_negatives'] + truths[2]

          print "\n\n -- RESULTS -- \n\n"
          # We've now built up the results dictionary for our K scores for all the documents.
          # Sum all the scores up, calculate their P/R/F and print it out
          for i in range(1, worker_k + 1):
            tp = len(results[i]['true_positives'])
            fp = len(results[i]['false_positives'])
            fn = len(results[i]['false_negatives'])

            precision, recall, f = self.determine_f( tp, fp, fn )

            # Make the summary files for the different K values
            with open('mturk_'+ disease_type +'_k'+ str(i) +'_summary.csv', 'wb') as csvfile:
              writer = csv.writer(csvfile, delimiter=',')
              writer.writerow(["TP", "FP", "FN", "precision", "recall", "F", "consistency"])

              arr = [tp, fp, fn, precision, recall, f, 100*2*float(tp)/(tp+tp+fp+fn)]
              writer.writerow(arr)
              print arr

            # Make the count files for the differnt K values
            for x in ['true_positives', 'false_positives', 'false_negatives']:
              with open('mturk_'+ disease_type +'_k'+ str(i) +'_'+ x +'.csv', 'wb') as csvfile:
                writer = csv.writer(csvfile, delimiter=',')
                writer.writerow(["text", "start", "count"])
                for term in results[i][x]:
                  arr = [term[0], term[1], ""]
                  writer.writerow(arr)



    def util_ncbo_specturm(self, documents):
        # ncbos = User.objects.filter(userprofile__ncbo = True).all()
        ncbos = Ncbo.objects.filter(score__lt = 15).all()
        gm_user = User.objects.get(username__exact = 'goldenmaster')

        with open('ncbo_spectrum.csv', 'wb') as csvfile:
          writer = csv.writer(csvfile, delimiter=',')
          writer.writerow(["Score", "Min Term Size", "TP", "FP", "FN", "P", "R", "F"])

          for ncbo in ncbos:
            results = []
            for document in documents:
              for section in document.section_set.all():
                # Collect the list of Annotations models for the Golden Master and NCBO Annotator to use throughout
                gm_annotations = list( Annotation.objects.filter(view__section = section, view__user = gm_user).values('text', 'start') )
                ncbo_annotations = list( Annotation.objects.filter(view__section = section, view__user = ncbo.user).values('text', 'start') )

                #
                # LOGGING
                #
                # print "\nGM Anns: "
                # for gm_ann in gm_annotations:
                #   print gm_ann.text + " :: " + str(gm_ann.start)

                # print "\nNCBO Anns: "
                # for ncbo_ann in ncbo_annotations:
                #   print ncbo_ann.text + " :: " + str(ncbo_ann.start)

                # print "\n - - - - - - \n"

                ncbo_score = self.calc_score(ncbo_annotations, gm_annotations)
                ncbo_score = ( len(ncbo_score[0]), len(ncbo_score[1]), len(ncbo_score[2]) )
                results.append( ncbo_score )


            results = map(sum,zip(*results))
            score = self.determine_f( results[0], results[1], results[2] )
            arr = [ncbo.score, ncbo.min_term_size, results[0], results[1], results[2], score[0], score[1], score[2]]
            writer.writerow(arr)
            print arr


    def match_exact(self, gm_ann, user_anns):
        '''
          Exact or Coextensive match finding for annotations. Works off start of annotation and cleaned length both being equal

          Returns True is any of the user annotations are equal to this GM Annotation

        '''
        gm_len = len(gm_ann['text'])
        for user_ann in user_anns:
          if gm_ann['start'] == user_ann['start'] and gm_len == len(user_ann['text']): return True
        return False


    def dict_to_tuple(self, dic):
        return (dic['text'], int(dic['start']))


    def calc_truth_table(self, annotations_a, annotations_b):
      '''
        This calculates the comparsion overlap between two arrays of dictionary terms

        It considers both the precision p and the recall r of the test to compute the score:
        p is the number of correct results divided by the number of all returned results
        r is the number of correct results divided by the number of results that should have been returned.
        The F1 score can be interpreted as a weighted average of the precision and recall, where an F1 score reaches its best value at 1 and worst score at 0.

       tp  fp
       fn  *tn

      '''
      true_positives = [gm_ann for gm_ann in annotations_b if self.match_exact(gm_ann, annotations_a)]

      # In order to make our comparisons we need to do the conversion to a list of tuples
      true_positives = [self.dict_to_tuple(item) for item in true_positives]
      annotations_a  = [self.dict_to_tuple(item) for item in annotations_a]
      annotations_b  = [self.dict_to_tuple(item) for item in annotations_b]

      # Annotations the user submitted that were wrong (the User set without their True Positives)
      false_positives = set(annotations_a) - set(true_positives)

      # Annotations the user missed (the GM set without their True Positives)
      false_negatives = set(annotations_b) - set(true_positives)

      # Add the False Positives and False Negatives to a global array to keep track which
      # annotations are commonly incorrect
      # for tp in true_positives:  self.error_aggreements['true_positives' ].append( tp[1][1] )
      # for fp in false_positives: self.error_aggreements['false_positives'].append( fp[1][1] )
      # for fn in false_negatives: self.error_aggreements['false_negatives'].append( fn[1][1] )

      return ( list(true_positives), list(false_positives), list(false_negatives) )


    def determine_f(self, true_positive, false_positive, false_negative):
        if true_positive + false_positive is 0:
          return (0,0,0)

        precision = true_positive / float(true_positive + false_positive)
        recall = true_positive / float(true_positive + false_negative)

        if precision + recall > 0.0:
          f = ( 2 * precision * recall ) / ( precision + recall )
          return (precision, recall, f)
        else:
          return (0,0,0)


    def compare_turk_to_gold():
      pass


    def get_user_data(self, username, source):
      views = View.objects.filter(user__username = username, section__document__source = source).all()
      return views


    def get_experiment_data(self, experiment):
      '''
        Returns all the Annotations for an experiment
        (TODO) Confirm they source documents for the correct experiment GM set
      '''
      if experiment is 1:
          # filter( Annotation.user.has(mturk=1) ).\
          # filter( Annotation.created < '2013-11-6' ).\
          annotations = Annotation.objects.filter(experiment = experiment).all()
      elif experiment is 2:
            # filter( Annotation.user.has(mturk=1) ).\
            # filter_by( experiment = experiment ).\
          annotations = Annotation.objects.filter(experiment = experiment).all()
      else:
          # annotations = Annotation.objects.filter(view__user__userprofile__mturk = True).all()
          annotations = Annotation.objects.filter(experiment = 3).all()
          views = View.objects.filter(user__userprofile__mturk = True).all()

      return views, annotations


