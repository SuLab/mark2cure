from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from mark2cure.document.models import Document, Section, View, Annotation

class Command(BaseCommand):
    args = '<experiment_run_id>'
    help = 'Run and calculate metrics for the experiment'

    def handle(self, *args, **options):
        for experiment in args:
            self.stdout.write('-- Running analytics on Experiment Run %s --' % experiment)

            # views, hits = self.get_experiment_data(experiment)
            gm_views = self.get_user_data('goldenmaster', 'NCBI_corpus_development')
            print gm_views
            print len(gm_views)
            counter = 0
            for gmv in gm_views:
              counter += len(gmv.annotation_set.all())

            print counter

            gm_anns = Annotation.objects.filter(user_agent = 'goldenmaster').all()
            print gm_anns
            print len(gm_anns)

            # for gmv in gm_views:
            #   gmv.delete()


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

