# -*- coding: utf-8 -*-
"""
    mark2cure.manage.analysis
    ~~~~~~~~~~~~~~~~~~~~~

    results/analysis commands
"""

from flask.ext.script import Command, prompt, prompt_pass
from mark2cure.settings import *
from ..core import db
from ..models import User, Document, Annotation, View

import requests, re

def user_vs_gold(user, doc):
    # This is a document that requires validation
    user_annotations = db.session.query(Annotation).filter_by(document = doc).filter_by(user = user).all()
    gold_annotations = db.session.query(Annotation).filter_by(document = doc).filter_by(user_id = 2).all()

    user_annotations = [ann.compare_view() for ann in user_annotations]
    gold_annotations = [ann.compare_view() for ann in gold_annotations]

    user_matches = len([ann for ann in user_annotations if ann in gold_annotations])

    return user_matches, ( len(user_annotations) - user_matches), len(user_annotations), len(gold_annotations)

class Compare(Command):
    "F Score"

    def run(self):
      document = db.session.query(Document).get(28)
      annotator_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user_id = 1).all()
      gold_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user_id = 2).all()

      print [ann.compare_view() for ann in annotator_annotations]
      print
      print [ann.compare_view() for ann in gold_annotations]


