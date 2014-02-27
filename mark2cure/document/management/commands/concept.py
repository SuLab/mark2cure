from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
# from django.db.models import Count
from django.contrib.auth.models import User
from django.utils.text import slugify

from mark2cure.document.models import *
from mark2cure.document.utils import create_from_pubmed_id

import os, os.path, csv
import warnings

class Command(BaseCommand):
    args = '<command>'
    help = 'Tools to import and manage concept and concept relationships'


    def handle(self, *args, **options):
        command = args[0]

        self.stdout.write('-- Running Concept Tools ({0}) --'.format(command))

        if command == "import":
          self.import_tsv()

        elif command == "average_time_per_hit":
          pass

        else:
          pass


    def import_tsv(self):
        with open('assets/datasets/concepts/chordoma.txt','r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader, None)  # skip the headers
            for pid, sid, pnumber, pmid, predicate, s_cui, s_name, s_type, s_novel, o_cui, o_name, o_type, o_novel in reader:
              print pid, sid, pnumber, pmid, predicate, s_cui, s_name, s_type, s_novel, o_cui, o_name, o_type, o_novel

              doc = create_from_pubmed_id(pmid)

              if doc:
                overview = doc.section_set.filter(kind = 'o').first()

                subject_concept, sc = Concept.objects.get_or_create(concept_id = s_cui, preferred_name = s_name)
                object_concept, oc = Concept.objects.get_or_create(concept_id = o_cui, preferred_name = o_name)
                relationship_type, rc = RelationshipType.objects.get_or_create(full_name = predicate, type = slugify(unicode(predicate)) )

                user, uc = User.objects.get_or_create(username="semmed")
                view, vc = View.objects.get_or_create(section = overview, user = user)
                ann,ac = Annotation.objects.get_or_create(view = view, kind = 'r')

                ConceptRelationship.objects.get_or_create(
                    concept = subject_concept,
                    relationship = relationship_type,
                    target = object_concept,
                    annotation = ann)

