# -*- coding: utf-8 -*-
"""
    mark2cure.manage.aws
    ~~~~~~~~~~~~~~~~~~~~~

    turker commands
"""

from flask.ext.script import Command, prompt, prompt_pass
from mark2cure.settings import *
from boto.mturk.connection import MTurkConnection, ExternalQuestion
from boto.mturk.question import QuestionContent,Question,QuestionForm,Overview,AnswerSpecification,SelectionAnswer,FormattedContent,FreeTextAnswer

import requests, datetime

class Turk(Command):
  "Add a mturk to a document"

  # Utility Methods
  def make_connection(self):
      return MTurkConnection( aws_access_key_id = AWS_ACCESS_ID,
                              aws_secret_access_key = AWS_SECRET_KEY,
                              host = AWS_HOST)

  def delete_all(self):
      mtc = self.make_connection()
      for hit in mtc.get_all_hits():
          mtc.disable_hit( hit.HITId )

  def ban_user(self, worker_id, reason="Did not meet our standards."):
      mtc = self.make_connection()
      mtc.block_worker(worker_id, reason)

  def external_question(self, doc_id):
      # Would be cool to pull the length of the document to know the correct size of the window to show
      return ExternalQuestion("https://mark2cure.org/mturk/#/document/"+str(doc_id), 800)

  def make_qualification_test(self):
      mtc = self.make_connection()
      name = 'Disease Recognitiion'
      description = 'Select the correct answer to move on'
      status = 'Active'
      keywords = 'science, annotation, disease'
      test_duration = 60 * 2

      question_form = QuestionForm()

      overview = Overview()
      overview.append_field('Title', 'Select the option which contains only diseases from the follow paragraph.')

      qc1 = QuestionContent()
      qc1.append_field('Title','Colorectal cancer occurs when tumors form in the lining of the large intestine. The risk of developing colorectal cancer rises after age 50. You\'re also more likely to get it if you have colorectal polyps, a family history of colorectal cancer, ulcerative colitis or Crohn\'s disease, eat a diet high in fat, or smoke.')

      fta2 = FreeTextAnswer()

      q1 = Question(identifier = 'ann_selection',
                    content = qc1,
                    answer_spec = AnswerSpecification(fta2),
                    is_required = True)

      question_form.append(overview)
      question_form.append(q1)

      qual_test = mtc.create_qualification_type(
                              name,
                              description,
                              status, keywords=None,
                              # Leave None to only allow once
                              retry_delay = None,
                              test =  question_form,
                              answer_key = None,
                              answer_key_xml = None,
                              test_duration = test_duration,
                              auto_granted = False)
      return qual_test

  def make_qualification_score(self):
      create_qualification_type(
        name = 'Mark2Cure GM Score',
        description = '',
        status = 'Active',
        auto_granted = True,
        auto_granted_value = 2)

  # Actionable methods
  def hit_for_document(self, doc_id, max_assignments = 5, reward = 0.04, minutes = 8, title="Annotate scientific articles"):
      mtc = self.make_connection()
      description = ('Visit a website and highlight diseases that are present in a paragraph.')
      keywords = 'science, annotation, disease'

      hit = mtc.create_hit(
          hit_type = None,
          question = self.external_question(doc_id),
          hit_layout = None,
          lifetime = datetime.timedelta(7),
          max_assignments = max_assignments,
          title = title,
          description = description,
          keywords = keywords,
          reward = reward,
          duration = datetime.timedelta(minutes = minutes),
          approval_delay = None,
          annotation = None,
          questions = None,
          qualifications = self.make_qualification_test(),
          layout_params = None,
          response_groups = None
          )
      return hit

  def run(self):
    documents = db.session.query(Document).all()
    for doc in documents:
        self.hit_for_document(doc.id)

    # Click Done:
    #   Submits their Annotations to our server
    #   Hit GM API to get pass / result value
    #   Update their score: http://docs.aws.amazon.com/AWSMechTurk/latest/AWSMturkAPI/ApiReference_UpdateQualificationScoreOperation.html
