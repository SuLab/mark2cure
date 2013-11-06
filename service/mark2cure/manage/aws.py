# -*- coding: utf-8 -*-
"""
    mark2cure.manage.aws
    ~~~~~~~~~~~~~~~~~~~~~

    turker commands
"""

from flask.ext.script import Command, prompt, prompt_pass
from mark2cure.settings import *

from boto.mturk.connection import *
from boto.mturk.question import *
from boto.mturk.qualification import *
from ..models import *
from ..core import db
import requests, datetime, random

class Turk(Command):
  "Add a mturk to a document"

  def __init__(self):
    self.mtc = MTurkConnection( aws_access_key_id = AWS_ACCESS_ID,
                              aws_secret_access_key = AWS_SECRET_KEY,
                              host = AWS_HOST)

  # Utility Methods
  def disable_all(self):
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

      # Give a short overview and GIF
      overview = Overview()
      overview.append_field('Title', 'Instructions')
      overview.append(FormattedContent( '<p>You will be presented with paragraphs from the biomedical literature which we believe may help resolve some important medically related questions. The task is to highlight words or phrases in that text which are diseases. Terms like breast cancer, juvenile diabetes, and others can be very difficult for a computer to recognize and we think people can do a better job!</p>'
                                        '<h2>Steps</h2>'
                                        '<ol>'
                                          '<li>'
                                            '<h3>Click to highlight as many <strong>diseases</strong> as you can.</h3>'
                                            '<img alt="Example of highlight behavior" src="https://mark2cure.org/img/instructions/1.gif" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                            '<h3>Highlight each time a disease is mentioned.</h3>'
                                            '<img alt="Example of highlighting each disease mention behavior" src="https://mark2cure.org/img/instructions/2.gif" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                            '<h3>Highlight diseases when mentioned as the <strong>full name</strong>, <strong>abbreviation</strong> or <strong>acronym</strong>.</h3>'
                                            '<img alt="Example of highlighting full name, abbreviations, and acronyms" src="https://mark2cure.org/img/instructions/3.gif" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                            '<h3>For diseases that span multiple words, <strong>click-and-drag</strong> to highlight all words together rather than clicking the words individually.</h3>'
                                            '<img alt="Example of highlighting using click-and-drag behavior" src="https://mark2cure.org/img/instructions/4.gif" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                        '</ol>'
                                        '<p>Use <strong>Google</strong>, <strong>Wikipedia</strong> or <strong>other resources</strong> to differentiate between diseases and other related words (like symptoms and genes).</p>'))

      # Define question content
      qc = QuestionContent()
      qc.append_field('Title', 'Select which of the options contains the most disease terms from the following sentence.')
      qc.append_field('Text', 'Colorectal cancer occurs when tumors form in the lining of the large intestine. The risk of developing colorectal cancer rises after age 50. You\'re also more likely to get it if you have colorectal polyps, a family history of colorectal cancer, ulcerative colitis or Crohn\'s disease, eat a diet high in fat, or smoke.')

      # Make question choices
      s1 = ("Colorectal cancer, Colorectal polyps, Crohn\'s disease", "A")
      s2 = ("Colorectal cancer, Polyps, Smoke, Lining", "B")
      s3 = ("Colorectal cancer, Ulcerative colitis, Crohn\'s disease", "C")
      s4 = ("Cancer, Polyps, Ulcerative colitis", "D")

      choices = SelectionAnswer(
          style='radiobutton',
          selections=[s1,s2,s3,s4])

      # Define question
      q = Question(identifier = 'ann_selection',
                    content = qc,
                    answer_spec = AnswerSpecification(choices),
                    is_required = True)

      # Add the content to the questionform
      question_form.append(overview)
      question_form.append(q)

      # Define evaluation mechanism
      answer_logic = '''<AnswerKey xmlns="http://mechanicalturk.amazonaws.com/AWSMechanicalTurkDataSchemas/2005-10-01/AnswerKey.xsd">
                          <Question>
                          <QuestionIdentifier>ann_selection</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>A</SelectionIdentifier>
                              <AnswerScore>0</AnswerScore>
                            </AnswerOption>
                            <AnswerOption>
                              <SelectionIdentifier>B</SelectionIdentifier>
                              <AnswerScore>0</AnswerScore>
                            </AnswerOption>
                            <AnswerOption>
                              <SelectionIdentifier>C</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                            <AnswerOption>
                              <SelectionIdentifier>D</SelectionIdentifier>
                              <AnswerScore>0</AnswerScore>
                            </AnswerOption>
                          </Question>
                        </AnswerKey>'''

      # qual_test = self.mtc.update_qualification_type(AWS_QUAL_TEST_2,
      qual_test = self.mtc.create_qualification_type(
        name = 'Disease recognition and UI Steps',
        description = 'Instructions and main steps to correctly use the highlighting behavior. Simple multiple-choice question to determine if the Worker understands the problem and has disease annotation ability.',
        status = 'Active',
        test = question_form,
        answer_key = answer_logic,
        # retry_delay = None,
        test_duration = 5 * 60)

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
    # score = self.mtc.update_qualification_type(AWS_QUAL_GM_SCORE,
    score = self.mtc.create_qualification_type(
        name = 'Golden master performance score',
        description = 'The score value which repersents how well a user annotates docuates with gm_validation=True',
        status = 'Active',
        auto_granted = True,
        auto_granted_value = 1)
    return score

  # Actionable methods
  def hit_for_document(self, doc_id, max_assignments = 5, reward = 0.06, minutes = 4, title="Highlight diseases in paragraph"):
      description = ('Highlight by clicking or dragging over multiple words in the following paragraph that are diseases. Do *not* select symptoms, conditions or any other non-disease term. When available, highlight multi-word disease together by click and dragging. When you accept the HIT, you will be allowed to start highlighting and a submit button will appear.')
      keywords = 'science, annotation, disease, text, highlight, annotation, medicine, term recognition'

      qualifications = Qualifications()
      # Add the simple test
      qualifications.add( Requirement(AWS_QUAL_TEST, "EqualTo", 1) )
      # Add the score
      # qualifications.add( Requirement(AWS_QUAL_GM_SCORE, "NotEqualTo", 0) )

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
    documents = db.session.query(Document).filter_by(source = 'NCBI_corpus_development').all()
    doc_ids = [doc.id for doc in documents]
    random.shuffle(doc_ids)
    for doc_id in doc_ids:
      print doc_id
      self.hit_for_document( doc_id )
