# -*- coding: utf-8 -*-
"""
    overholt.manage.setup
    ~~~~~~~~~~~~~~~~~~~~~

    content prep / setting up document commands
"""

from flask.ext.script import Command, prompt, prompt_pass
from ..models import *
from mark2cure.settings import *

import requests

class Heatmap(Command):
    "Calculates the heatmap array"

    def run(self):
      for document in Document.query.all():
        print "Running document {0}".format( document.document_id )
        # Select all the annotations for this document
        annotations = db.session.query(Annotation).filter_by(document = document)
        step = 0
        length = 0
        pop_arr = []

        # Iterate over each word in the text
        for idx, word in enumerate(document.text.split()):
          length = len(word)
          step = step + length + 1;
          wordstart = step - length - 1,

          # Count which words have annotations in the db for that doc, across all users
          # Number of annotations that encapsulate that word
          count = 0
          for ann in annotations:
            if wordstart[0] >= ann.start and wordstart[0] <= (ann.start+ann.length):
              count = count + 1
          pop_arr.append(count);

        document.cache = ', '.join([str(x) for x in pop_arr])

        # Save every document instead of once incase some doc crashes
        db.session.commit()
      print "Complete"

class Annotate(Command):
    "Populate the documents with the NCBI Annotator"

    def run(self):
      # Get the annotator bot account to make the db entries
      user = db.session.query(User).get(1)

      # NCBO Annotator http request
      payload = { 'apikey'                    : NCBO_API_KEY,
                  'stopWords'                 : STOP_WORDS,
                  'minTermSize'               : 5,
                  'withSynonyms'              : 'true',
                  # 1009 = Human disease ontology
                  'ontologiesToKeepInResult'  : '1009',
                  'isVirtualOntologyId'       : 'true',
                  # T047 -- Disease or Syndrome
                  # 'semanticTypes'             : 'T047',
                  'textToAnnotate'            : '',
                }

      # Go over all the documents
      documents = Document.query.all()
      for document in documents:
        # If the current document doesn't have at least 4 annotations from the bot, try to get more...
        if db.session.query(Annotation).filter_by(document = document).filter_by(user = user).count() <= 3:
          print "Running document {0}".format( document.document_id )
          payload['textToAnnotate'] = document.text
          r = requests.post("http://rest.bioontology.org/obs/annotator", data=payload)

          if r.ok:
            try:
              root = ET.fromstring( r.text.decode('utf-8') )
              # Make array of all relevant NCBO results
              for ann in root.iter('annotationBean'):
                ctx = ann.find('context')
                start = int(ctx.find('from').text)-1
                stop = int(ctx.find('to').text)
                annotation = document.text[start:stop]
                concept_url =  ctx.find('term').find('concept').find('fullId').text

                concept = db.session.query(Concept).filter_by(concept_id = concept_url).first()
                if concept is None:
                  concept = Concept(concept_url)
                  db.session.add(concept)
                  db.session.commit()

                # Find those results in the original abstract we had
                ann = Annotation( 0,
                                  'disease',
                                  annotation,
                                  start,
                                  len(annotation),
                                  stop,
                                  user,
                                  document,
                                  'annotator_bot_1.0',
                                  '',
                                  concept
                                );
                db.session.add(ann)

              # Save every document instead of once incase some doc crashes
              db.session.commit()
            except Exception:
              pass
      print ":)"

