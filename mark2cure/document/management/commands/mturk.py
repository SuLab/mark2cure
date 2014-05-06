from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.contrib.auth.models import User

from mark2cure.document.models import Document, Section, View, Annotation
from mark2cure.account.models import Ncbo

from mark2cure.common.aws import Turk

class Command(BaseCommand):
    args = '<experiment_run_id>'
    help = 'Command for posting and controlling mturk activity'

    def handle(self, *args, **options):
        if len(args) < 1: raise Exception('Analysis needs 2 parameters <experiement_id, command>')

        command = args[0]
        hit_count = 0
        if len(args) > 1:
            hit_count = args[1]

        self.stdout.write('-- Running MTurk Commands ({0}) on Documents {1} --'.format(command, document_set))

        turk = Turk()

        if command == "create_qual":
            print turk.make_qualification_test()


        elif command == "disable_all":
            turk.disable_all()


        elif command == "create_hits":
            for idx in range(hit_count):
                turk.hit_for_document(max_assignments = 1, minutes = 10)
                print idx


        else:
          pass


