from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from mark2cure.document.models import Document, Section, View, Annotation

import os, os.path, csv


class Command(BaseCommand):
    args = '<experiment_run_id>'
    help = 'Run and calculate metrics for the experiment'

    def handle(self, *args, **options):
        for experiment in args:
            self.stdout.write('-- Running analytics on Experiment Run %s --' % experiment)

            # views, hits = self.get_experiment_data(experiment)
            # gm_views = self.get_user_data('goldenmaster', 'NCBI_corpus_development')
            # gm_anns = Annotation.objects.filter(user_agent = 'goldenmaster').all()

            # gm_counter = 0
            # for v in gm_views:
            #   gm_counter += len( v.annotation_set.all() )


            # print gm_counter
            # print len(gm_anns)


            types = ["disease:modifier", "disease:class", "disease:specific", "disease:composite"]
            for disease_type in types:
              gm_anns = Annotation.objects.filter(user_agent = 'goldenmaster', type=disease_type).all()
              print self.calc_score(gm_anns, gm_anns)

              # print "\n\n"+ disease_type +"\n\n"
              # for gm_ann in gm_anns[:10]:
                # print gm_ann.start



            # base_path = os.path.join(settings.PROJECT_PATH, 'results')
            # os.chdir(base_path)

            # Annotator vs. Gold

        #     max_k = 0

        #     for ann in mturk_anns:
        #

        #       for k in range(1, 5+1):

        #

        #     types = ["disease:modifier", "disease:class", "disease:specific", "disease:composite"]
        #     for disease_type in types:
        #       report = compareAnnosCorpusLevel(gold_annos, test_annos);




        # if(include_doc_level) {
        #   List<Report> reports = compareAnnosDocumentLevel(gold_annos, test_annos);
        #   report.writeReportList(reports, outprefix+"_doc_level.txt");
        # }


        #     with open('.csv', 'wb') as csvfile:
        #       writer = csv.writer(csvfile, delimiter=',')
        #       writer.writerow(['foo', 'bar', 'tall'])


        #     # print len(views)
        #     # print len(gm_views)
        # pass


    def match_exact(self, gm_ann, user_anns):
        '''
          Exact or Coextensive match finding for annotations. Works off start of annotation and cleaned length both being equal

          Returns True is any of the user annotations are equal to this GM Annotation

        '''
        gm_len = len(gm_ann.text)
        for user_ann in user_anns:
          if gm_ann.start == user_ann.start and gm_len == len(user_ann.text): return True
        return False



    def calc_score(self, annotations_a, annotations_b):
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
      true_positives = [item.simple() for item in true_positives]
      annotations_a  = [item.simple() for item in annotations_a]
      annotations_b  = [item.simple() for item in annotations_b]

      # Annotations the user submitted that were wrong (the User set without their True Positives)
      false_positives = set(annotations_a) - set(true_positives)

      # Annotations the user missed (the GM set without their True Positives)
      false_negatives = set(annotations_b) - set(true_positives)

      # Add the False Positives and False Negatives to a global array to keep track which
      # annotations are commonly incorrect
      # for tp in true_positives:  self.error_aggreements['true_positives' ].append( tp[1][1] )
      # for fp in false_positives: self.error_aggreements['false_positives'].append( fp[1][1] )
      # for fn in false_negatives: self.error_aggreements['false_negatives'].append( fn[1][1] )

      return ( len(true_positives), len(false_positives), len(false_negatives) )


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

