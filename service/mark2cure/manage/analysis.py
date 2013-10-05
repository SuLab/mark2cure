# -*- coding: utf-8 -*-
"""
    mark2cure.manage.analysis
    ~~~~~~~~~~~~~~~~~~~~~

    results/analysis commands
"""

from flask.ext.script import Command, prompt, prompt_pass
from mark2cure.settings import *

import requests, re


class Compare(Command):
    "F Score"

    def run(self):
      document = db.session.query(Document).get(28)
      annotator_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user_id = 1).all()
      gold_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user_id = 2).all()

      print [ann.compare_view() for ann in annotator_annotations]
      print
      print [ann.compare_view() for ann in gold_annotations]


