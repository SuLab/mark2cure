# -*- coding: utf-8 -*-
"""
    mark2cure.manage.aws
    ~~~~~~~~~~~~~~~~~~~~~

    turker commands
"""

from flask.ext.script import Command, prompt, prompt_pass

import settings, requests, re, csv, datetime

from boto.mturk.connection import MTurkConnection, ExternalQuestion



class Turk(Command):
  "Add a mturk to a document"

  def run(self):
    mtc = MTurkConnection(aws_access_key_id = settings.AWS_ACCESS_ID,
                          aws_secret_access_key = settings.AWS_SECRET_KEY,
                          host = settings.HOST)
    title = 'Help annotate scientific articles'
    description = ('Visit a website and highlight diseases that are present in a paragraph.')
    keywords = 'science, annotation, disease'

    q = ExternalQuestion("http://beta.mark2cure.org/#/22", 800)
    hit = mtc.create_hit(question = q,
               max_assignments = 1,
               title = title,
               description = description,
               keywords = keywords,
               duration = 60*5,
               reward = 0.01)

    # http://beta.mark2cure.org/#/22?
    #    assignmentId=2FKW6NG1FS3N2YP7T6SF2SZL3JWKY2
    #    &hitId=22DWJ5OPB0YVKAEUHN27HP9UE9K5XV
    #    &workerId=A296YJ2WQNOSKY
    #    &turkSubmitTo=https%3A%2F%2Fworkersandbox.mturk.com
    return hit

