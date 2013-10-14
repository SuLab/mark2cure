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
from boto.mturk.qualification import AdultRequirement, LocaleRequirement, NumberHitsApprovedRequirement, PercentAssignmentsAbandonedRequirement, PercentAssignmentsApprovedRequirement, PercentAssignmentsRejectedRequirement, PercentAssignmentsReturnedRequirement, PercentAssignmentsSubmittedRequirement, Qualifications, Requirement

import requests, datetime

class Turk(Command):
  "Add a mturk to a document"

  def __init__(self):
    self.mtc = MTurkConnection( aws_access_key_id = AWS_ACCESS_ID,
                              aws_secret_access_key = AWS_SECRET_KEY,
                              host = AWS_HOST)

  # Utility Methods
  def delete_all(self):
      for hit in self.mtc.get_all_hits():
          self.mtc.disable_hit( hit.HITId )

  def ban_user(self, worker_id, reason="Did not meet our standards."):
      self.mtc.block_worker(worker_id, reason)

  def external_question(self, doc_id):
      # Would be cool to pull the length of the document to know the correct size of the window to show
      return ExternalQuestion("https://mark2cure.org/mturk/#/document/"+str(doc_id), 800)

  def make_qualification_test(self):
      '''
        This qualification type only gets created once and is what is shown to
        new workers before answering any questions

      '''
      # Define question form
      question_form = QuestionForm()

      # Define question content
      qc = QuestionContent()
      qc.append_field('Title', 'Select which of the follow options contains the most disease terms')
      qc.append_field('Text', 'Colorectal cancer occurs when tumors form in the lining of the large intestine. The risk of developing colorectal cancer rises after age 50. You\'re also more likely to get it if you have colorectal polyps, a family history of colorectal cancer, ulcerative colitis or Crohn\'s disease, eat a diet high in fat, or smoke.')

      # Make question choices
      s1 = '<Selection><SelectionIdentifier>s1</SelectionIdentifier>AAAAAA</Selection>'
      s2 = '<Selection><SelectionIdentifier>s2</SelectionIdentifier>BBBBBB</Selection>'
      s3 = '<Selection><SelectionIdentifier>s3</SelectionIdentifier>CCCCCC</Selection>'
      s4 = '<Selection><SelectionIdentifier>s4</SelectionIdentifier>DDDDDD</Selection>'
      # s1 = "Colorectal cancer, Crohn\'s disease, colorectal polyps"
      # s2 = "Polyps, smoke, lining, Colorectal cancer"
      # s3 = "Colorectal cancer, ulcerative colitis, Crohn\'s disease"
      # s4 = "Ulcerative colitis, polyps, cancer"

      choices = SelectionAnswer(
          style='multichooser',
          selections=[s1,s2,s3,s4])

      # Define question
      q = Question(identifier = 'ann_selection',
                    content = qc,
                    answer_spec = AnswerSpecification(choices),
                    is_required = True)

      question_form.append(q)

      # Define evaluation mechanism
      answer_logic = '''<Question>
                        <QuestionIdentifier>ann_selection</QuestionIdentifier>
                        <AnswerOption>
                          <SelectionIdentifier>s3</SelectionIdentifier>
                          <AnswerScore>1</AnswerScore>
                        </AnswerOption>
                      </Question>'''

      qual_test = self.mtc.create_qualification_type(
        name = 'mark2cure_simple_recognitiion_1',
        description = 'Simple multiple-choice form to determine if the Worker understands the problem and has basic disease annotation ability',
        status = 'Active',
        # Leave None to only allow once
        retry_delay = None,
        test =  question_form,
        answer_key_xml = answer_logic,
        test_duration = 60 * 2,
        auto_granted = False)
      return qual_test

  def make_qualification_score(self):
    '''
      This qualification type is what allows us to use Golden Master score values checks on defined
      documents used to score.

      0 -- Booted from HIT submission
      1 -- Required to play
      2 -- Starting value
      3
    '''
    score = self.mtc.create_qualification_type(
        name = 'mark2cure_gm_score',
        description = 'The score value which repersents how well a user annotates docuates with gm_validation=True',
        status = 'Active',
        auto_granted = True,
        auto_granted_value = 2)
    return score

  def score_qualification(self, qualification_type_id, worker_id):
    pass

  def intro_test_qualification(self):
    pass

  # Actionable methods
  def hit_for_document(self, doc_id, max_assignments = 5, reward = 0.04, minutes = 8, title="Annotate scientific articles"):
      description = ('Visit a website and highlight diseases that are present in a paragraph.')
      keywords = 'science, annotation, disease'

      qualifications = Qualifications()
      qualifications.add(  self.intro_test_qualification() )

      hit = self.mtc.create_hit(
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
          qualifications = qualifications,
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
