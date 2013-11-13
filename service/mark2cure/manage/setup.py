# -*- coding: utf-8 -*-
"""
    overholt.manage.setup
    ~~~~~~~~~~~~~~~~~~~~~

    content prep / setting up document commands
"""

from flask.ext.script import Command, prompt, prompt_pass
from ..models import *
from ..core import db

from mark2cure.settings import *

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


