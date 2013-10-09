# -*- coding: utf-8 -*-
"""
    mark2cure.manage.aws
    ~~~~~~~~~~~~~~~~~~~~~

    turker commands
"""

from flask.ext.script import Command, prompt, prompt_pass
from mark2cure.settings import *
from boto.mturk.connection import MTurkConnection, ExternalQuestion

import requests


def make_connection():
  return MTurkConnection(aws_access_key_id = AWS_ACCESS_ID,
                          aws_secret_access_key = AWS_SECRET_KEY,
                          host = AWS_HOST)

def clear_all():
  mtc = make_connection()
  for hit in mtc.get_all_hits():
      mtc.disable_hit( hit.HITId )


class Turk(Command):
  "Add a mturk to a document"

  def run(self):
    mtc = self.make_connection()
    title = 'Help annotate scientific articles'
    description = ('Visit a website and highlight diseases that are present in a paragraph.')
    keywords = 'science, annotation, disease'

    q = ExternalQuestion("https://mark2cure.org/mturk/#/document/25", 600)
    hit = mtc.create_hit(question = q,
               max_assignments = 5,
               title = title,
               description = description,
               keywords = keywords,
               duration = 60*8,
               reward = 0.04)



    q = ExternalQuestion("http://beta.mark2cure.org/#/document/22", 800)

    # Click Done:
    #   Submits their Annotations to our server
    #   Hit GM API to get pass / result value
    #   Update their score: http://docs.aws.amazon.com/AWSMechTurk/latest/AWSMturkAPI/ApiReference_UpdateQualificationScoreOperation.html

    return hit

