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

def make_qualification_test():
  create_qualification_type(name,
                            description,
                            status, keywords=None,
                            retry_delay=None,
                            test=None,
                            answer_key=None,
                            answer_key_xml=None,
                            test_duration=None,
                            auto_granted=False,
                            auto_granted_value=1)


class Turk(Command):
  "Add a mturk to a document"

  def run(self):
    mtc = self.make_connection()
    title = 'Help annotate scientific articles'
    description = ('Visit a website and highlight diseases that are present in a paragraph.')
    keywords = 'science, annotation, disease'

    q = ExternalQuestion("https://mark2cure.org/mturk/#/document/16", 800)
    hit = mtc.create_hit(question = q,
               max_assignments = 5,
               title = title,
               description = description,
               keywords = keywords,
               duration = 60*8,
               reward = 0.04)


    hit = mtc.create_hit( hit_type=None, 
        question=None, 
        hit_layout=None, 
        lifetime=datetime.timedelta(7), 
        max_assignments=1, 
        title=None, 
        description=None, 
        keywords=None, 
        reward=None, 
        duration=datetime.timedelta(7), 
        approval_delay=None, 
        annotation=None, 
        questions=None, 
        qualifications=None, 
        layout_params=None, 
        response_groups=None)¶



    # Click Done:
    #   Submits their Annotations to our server
    #   Hit GM API to get pass / result value
    #   Update their score: http://docs.aws.amazon.com/AWSMechTurk/latest/AWSMturkAPI/ApiReference_UpdateQualificationScoreOperation.html

    return hit

