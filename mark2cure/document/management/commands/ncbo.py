from django.conf import settings
from django.core.management.base import BaseCommand, CommandError, NoArgsCommand
from django.contrib.auth.models import User

from mark2cure.document.models import Document, Section, View, Annotation, Concept
from mark2cure.account.models import Ncbo

import requests, re, csv, datetime
import xml.etree.ElementTree as ET


class Command(BaseCommand):
    help = 'Populates the documents with NCBO annotator results'

    def handle(self, *args, **options):
        # self.stdout.write('-- Populating documents with NCBO annotator results --')

        for minsize in range(3, 8):
          for score in range(5, 80, 10):
            name = 'ncbo_'+ str(minsize) +'_'+ str(score)

            u, uc = User.objects.get_or_create(username = name)
            if uc:
              u.set_password('')
              profile = u.profile
              profile.email_notify = False
              profile.ncbo = True
              profile.save()
              u.save()

            ncbo, created = Ncbo.objects.get_or_create(min_term_size=minsize, score=score, user=u)

        # Go over all the documents
        annotators = Ncbo.objects.filter(score = 5).all()
        documents = Document.objects.all()

        for document in documents:
          for annotator in annotators:
            # If the current document doesn't have at least 4 annotations from the bot, try to get more...
            current_anns = Annotation.objects.filter(view__section__document = document, view__user = annotator.user ).count()

            if current_anns <= 3:
              for section in document.section_set.all():
                view, vc = View.objects.get_or_create(section = section, user = annotator.user)

                # NCBO Annotator http request
                payload = { 'apikey'             : settings.NCBO_API_KEY,
                            'stopWords'                 : settings.STOP_WORDS,
                            'minTermSize'               : 5,
                            'withSynonyms'              : 'true',
                            # 1009 = Human disease ontology
                            'ontologiesToKeepInResult'  : '1009',
                            'isVirtualOntologyId'       : 'true',
                            # T047 -- Disease or Syndrome
                            # 'semanticTypes'             : 'T047',
                            'textToAnnotate'            : '',
                        }


                # Config settings to change per request
                payload['textToAnnotate'] = view.section.text
                payload['minTermSize'] = annotator.min_term_size

                r = requests.post("http://rest.bioontology.org/obs/annotator", data=payload)
                self.handle_request(r, annotator, view)

    def handle_request(self, request, annotator, view):
      if request.ok:
        try:
          root = ET.fromstring( request.text.decode('utf-8') )

          # Make array of all relevant NCBO results
          for ann in root.iter('annotationBean'):
            ctx = ann.find('context')
            score = int(ann.find('score').text)
            start = int(ctx.find('from').text)-1
            stop = int(ctx.find('to').text)

            annotation = view.section.text[start:stop]
            concept_url =  ctx.find('term').find('concept').find('fullId').text
            concept_preferred_name = ctx.find('term').find('concept').find('preferredName').text

            print "\n - - - - - - CONCEPT - - - - - - - - - \n"
            print concept_preferred_name

            if score >= annotator.score:
              concept, cc = Concept.objects.get_or_create(full_id = concept_url)
              concept.preferred_name = concept_preferred_name
              concept.save()

              # Find those results in the original abstract we had
              ann = Annotation()
              ann.kind = 'e'
              ann.type = 'disease'
              ann.text = annotation
              ann.start = start
              ann.user_agent = 'ncbo'
              ann.view = view
              ann.concept = concept
              ann.save()

        except Exception:
          pass
