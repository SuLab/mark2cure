from django.conf import settings
from datetime import datetime
from django.contrib.auth.models import User

from mark2cure.document.models import Document, Section, View, Annotation

from boto.mturk.connection import *
from boto.mturk.question import *
from boto.mturk.qualification import *

import requests, datetime, random, re, nltk

def get_mturk_account(worker_id):
    u, created = User.objects.get_or_create(username=worker_id)
    if created:
        u.set_password('')
        profile = u.profile
        profile.email_notify = False
        profile.mturk = True
        profile.save()
        u.save()

    return u



def get_timezone_offset(tz):
    offset = datetime.now(tz).strftime('%z')

    sign = 1
    if offset[0] == '-':
        sign = -1

    return sign * float(offset[1:]) / 100

'''
  Text Utility Functions

'''

def is_negative(text):
    return False


def text_has_link(text):
    pattern = re.compile("(https?://|a href)")
    match = pattern.findall(text, re.I)
    return True if len(match) else False


def text_has_any_words(text, words):
    for word in text.split():
        if word in words:
            return True

    return False


def text_has_any_string(text, strings):
    for string in strings:
        if string in text:
            return True

    return False


def text_has_any_pattern(text, patterns):
    combined = "(" + "|".join(patterns) + ")"
    if re.match(combined, text):
        return True
    return False


'''
  Language Utilities
'''

def tokenize(text):
    try:
        return nltk.word_tokenize(text)
    except Exception, e:
        print e
        return text.split()

def pos_tag(text):
    tokens = tokenize(text)
    try:
        return nltk.pos_tag(tokens)
    except Exception:
        tags = []
        for word in tokens:
            tags.append((word, '???'))

        return tags


def get_proper_nouns(text):
    pos_tags = pos_tag(text)
    return [word for word,pos in pos_tags if pos == 'NNP']


def language_clean(text):
    text = text.strip()
    # text.delete!("\"") if text.count("\"").odd? # derp, we have an unmatched quote? -- let's remove ALL the quotes!
    # text.gsub!(/\A\)[\.\s]/, "") # dumbass parens
    # text = text.strip.strip_html.squeeze(" ").gsub(/\.{4,}/,
    # "...").gsub(/(\r|\n|\t)/," ").gsub(/\?{2,}/, "?").gsub(/\!{2,}/, "!")
    return text


def clean_spaces(text):
    return ' '.join(text.split()).strip()


def is_question(text):
    return '?' in text.strip()[-3:]


def is_respondable(text):
    return True

'''
  AWS Stuff
  ---------

  Helper methods for creating HITs, Qualifications
'''

class Turk():

  def __init__(self):
    self.mtc = MTurkConnection( aws_access_key_id = settings.AWS_ACCESS_ID,
                                aws_secret_access_key = settings.AWS_SECRET_KEY,
                                host = settings.AWS_HOST)

  # Utility Methods
  def disable_all(self):
      for hit in self.mtc.get_all_hits():
          self.mtc.disable_hit( hit.HITId )

  def ban_user(self, worker_id, reason="Did not meet our standards."):
      self.mtc.block_worker(worker_id, reason)

  def external_question(self, doc_id):
      # Would be cool to pull the length of the document to know the correct size of the window to show
      return ExternalQuestion("https://mark2cure.org/document/"+str(doc_id)+"/", 800)

  def make_qualification_test(self):
      '''
        This qualification type only gets created once and is what is shown to
        new workers before answering any questions

      '''
      # Define question form
      question_form = QuestionForm()

      #
      # Instructions to train the Worker
      #
      overview = Overview()
      overview.append_field('Title', 'Instructions')
      overview.append(FormattedContent( '<p><strong>Task:</strong> You will be presented with text from the biomedical literature which we believe may help resolve some important medical related questions. The task is to highlight words and phrases in that text which are <u>diseases</u>, <u>disease groups</u>, or <u>symptoms</u> of diseases.  <u>This work will help advance research in cancer and many other diseases</u>!</p>'
                                        '<p><strong>Here are some examples of correctly highlighted text.  Please study these before attempting to take the qualification test.  Please also feel free to refer back to these examples if you are uncertain.</strong></p>'
                                        '<ul>'
                                          '<li>'
                                            '<h3>Rule #1: Highlight <u>all</u> diseases and disease abbreviations</h3>'
                                            '<img alt="Highlight all diseases and disease abbreviations" src="http://mark2cure.org/static/images/experiment/3/1.png" />'
                                            '<img alt="Highlight all diseases and disease abbreviations" src="http://mark2cure.org/static/images/experiment/3/2.png" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                            '<h3>Rule #2: Highlight the longest span of text specific to a disease</h3>'
                                            '<img alt="Highlight the longest span of text specific to a disease" src="http://mark2cure.org/static/images/experiment/3/3.png" />'
                                            '<img alt="Highlight the longest span of text specific to a disease" src="http://mark2cure.org/static/images/experiment/3/4.png" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                            '<h3>Rule #3: Highlight disease conjunctions as single, long spans.</h3>'
                                            '<img alt="Highlight disease conjunctions as single, long spans" src="http://mark2cure.org/static/images/experiment/3/5.png" />'
                                            '<img alt="Highlight disease conjunctions as single, long spans" src="http://mark2cure.org/static/images/experiment/3/6.png" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                          '<h3>Rule #4: Highlight symptoms - physical results of having a disease</h3>'
                                            '<img alt="Highlight symptoms - physical results of having a disease" src="http://mark2cure.org/static/images/experiment/3/7.png" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                            '<h3>Rule #5: Highlight <u>all</u> occurrences of disease terms</h3>'
                                            '<img alt="Highlight all occurrences of disease terms" src="http://mark2cure.org/static/images/experiment/3/8.png" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                            '<h3>Rule #6: Highlight <u>all</u> diseases, disease groups and key disease symptoms</h3>'
                                            '<img alt="Highlight all diseases, disease groups and key disease symptoms" src="http://mark2cure.org/static/images/experiment/3/9.png" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                        '</ul>'))

      #
      # Questions to ask the Worker
      #
      instructions = "Select all and only the terms that should be highlighted for each text segment (don't select terms that overlap with each other in the text):"
      paragraph1 = "Test Paragraph 1: \"Myotonic dystrophy ( DM ) is associated with a ( CTG ) n trinucleotide repeat expansion in the 3-untranslated region of a protein kinase-encoding gene , DMPK , which maps to chromosome 19q13 . 3 . \""
      paragraph2 = "Test Paragraph 2: \"Germline mutations in BRCA1 are responsible for most cases of inherited breast and ovarian cancer . However , the function of the BRCA1 protein has remained elusive . As a regulated secretory protein , BRCA1 appears to function by a mechanism not previously described for tumour suppressor gene products.\""
      paragraph3 = "Test Paragraph 3: \"We report about Dr . Kniest , who first described the condition in 1952 , and his patient , who , at the age of 50 years is severely handicapped with short stature , restricted joint mobility , and blindness but is mentally alert and leads an active life .  This is in accordance with molecular findings in other patients with Kniest dysplasia and...\""

      # # # # # # # # #
      #
      # Paragraph 1
      #
      # # # # # # # # #
      #abcdefghijklmnopqrstuvwxyz
      qc = QuestionContent()
      qc.append_field('Title', "Select the words that should be highlighted in the next three test paragraphs according to the rules above.")
      qc.append_field('Title', paragraph1)
      qc.append_field('Title', "")
      qc.append_field('Title', "Which of the following should be highlighted according to the instructions described above?")
      qc.append_field('Text', 'Myotonic')
      q1 = Question(identifier = 'term_selection_1', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t1"), ("False", "f1")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'dystrophy')
      q2 = Question(identifier = 'term_selection_2', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t2"), ("False", "f2")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'Myotonic dystrophy')
      q3 = Question(identifier = 'term_selection_3', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t3"), ("False", "f3")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'DM')
      q4 = Question(identifier = 'term_selection_4', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t4"), ("False", "f4")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'CTG')
      q5 = Question(identifier = 'term_selection_5', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t5"), ("False", "f5")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'trinucleotide repeat expansion')
      q6 = Question(identifier = 'term_selection_6', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t6"), ("False", "f6")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'kinase-encoding gene')
      q7 = Question(identifier = 'term_selection_7', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t7"), ("False", "f7")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'DMPK')
      q8 = Question(identifier = 'term_selection_8', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t8"), ("False", "f8")])))

      # # # # # # # # #
      #
      # Paragraph 2
      #
      # # # # # # # # #
      qc = QuestionContent()
      qc.append_field('Title', paragraph2)
      qc.append_field('Text', 'Germline mutations')
      q9 = Question(identifier = 'term_selection_9', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t9"), ("False", "f9")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'BRCA1')
      q10 = Question(identifier = 'term_selection_10', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t10"), ("False", "f10")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'breast')
      q11 = Question(identifier = 'term_selection_11', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t11"), ("False", "f11")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'ovarian cancer')
      q12 = Question(identifier = 'term_selection_12', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t12"), ("False", "f12")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'inherited breast and ovarian cancer')
      q13 = Question(identifier = 'term_selection_13', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t13"), ("False", "f13")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'cancer')
      q14 = Question(identifier = 'term_selection_14', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t14"), ("False", "f14")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'tumour')
      q15 = Question(identifier = 'term_selection_15', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t15"), ("False", "f15")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'tumour suppressor')
      q16 = Question(identifier = 'term_selection_16', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t16"), ("False", "f16")])))


      # # # # # # # # #
      #
      # Paragraph 3
      #
      # # # # # # # # #
      qc = QuestionContent()
      qc.append_field('Title', paragraph3)
      qc.append_field('Text', 'age of 50 years')
      q17 = Question(identifier = 'term_selection_17', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t17"), ("False", "f17")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'severely handicapped')
      q18 = Question(identifier = 'term_selection_18', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t18"), ("False", "f18")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'short')
      q19 = Question(identifier = 'term_selection_19', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t19"), ("False", "f19")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'short stature')
      q20 = Question(identifier = 'term_selection_20', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t20"), ("False", "f20")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'restricted joint mobility')
      q21 = Question(identifier = 'term_selection_21', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t21"), ("False", "f21")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'blindness')
      q22 = Question(identifier = 'term_selection_22', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t22"), ("False", "f22")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'mentally alert')
      q23 = Question(identifier = 'term_selection_23', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t23"), ("False", "f23")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'molecular findings')
      q24 = Question(identifier = 'term_selection_24', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t24"), ("False", "f24")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'Kniest dysplasia')
      q25 = Question(identifier = 'term_selection_25', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t25"), ("False", "f25")])))

      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'dysplasia')
      q26 = Question(identifier = 'term_selection_26', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t26"), ("False", "f26")])))

      # Add the content to the questionform
      question_form.append(overview)
      question_form.append(q1)
      question_form.append(q2)
      question_form.append(q3)
      question_form.append(q4)
      question_form.append(q5)
      question_form.append(q6)
      question_form.append(q7)
      question_form.append(q8)
      question_form.append(q9)
      question_form.append(q10)
      question_form.append(q11)
      question_form.append(q12)
      question_form.append(q13)
      question_form.append(q14)
      question_form.append(q15)
      question_form.append(q16)
      question_form.append(q17)
      question_form.append(q18)
      question_form.append(q19)
      question_form.append(q20)
      question_form.append(q21)
      question_form.append(q22)
      question_form.append(q23)
      question_form.append(q24)
      question_form.append(q25)
      question_form.append(q26)

      # Define evaluation mechanism
      answer_logic = '''<AnswerKey xmlns="http://mechanicalturk.amazonaws.com/AWSMechanicalTurkDataSchemas/2005-10-01/AnswerKey.xsd">

                          <Question>
                          <QuestionIdentifier>term_selection_1</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f1</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_2</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f2</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_3</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>t3</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_4</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>t4</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_5</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f5</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_6</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f6</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_7</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f7</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_8</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f8</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>





                          <Question>
                          <QuestionIdentifier>term_selection_9</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f9</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_10</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f10</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_11</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f11</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_12</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f12</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_13</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>t13</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_14</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f14</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_15</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>t15</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_16</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f16</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>



                          <Question>
                          <QuestionIdentifier>term_selection_17</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f17</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_18</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>t18</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_19</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f19</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_20</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>t20</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_21</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>t21</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_22</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>t22</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_23</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f23</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_24</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f24</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_25</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>t25</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_26</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>f26</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                        </AnswerKey>'''

      # qual_test = self.mtc.update_qualification_type(settings.AWS_QUAL_TEST_3,
      qual_test = self.mtc.create_qualification_type(
        name = 'Annotation Instructions & Qualification Questions',
        description = 'Detailed annotation instructions. Multiple-choice questions to assess concept understanding.',
        status = 'Active',
        test = question_form,
        answer_key = answer_logic,
        # retry_delay = 1,
        test_duration = 25 * 60)

      return qual_test

  # Actionable methods
  def hit_for_document(self, doc_id, max_assignments = 5, reward = 0.06, minutes = 4, title="Highlight diseases in paragraph"):
      description = ('You will be presented with text from the biomedical literature which we believe may help resolve some important medically related questions. The task is to highlight words and phrases in that text which are, or are highly related to, diseases.  This work will help advance research in cancer and many other diseases!')
      keywords = 'science, abstract, primary literature, annotation, disease, text, highlight, annotation, medicine, term recognition'

      qualifications = Qualifications()
      # Add the step instructions and basic test
      qualifications.add( Requirement(settings.AWS_QUAL_TEST_3, "GreaterThanOrEqualTo", 22) )

      hit = self.mtc.create_hit(
          hit_type = None,
          question = self.external_question(doc_id),
          hit_layout = None,
          lifetime = datetime.timedelta(14),
          max_assignments = max_assignments,
          title = title,
          description = description,
          keywords = keywords,
          reward = reward,
          duration = datetime.timedelta(minutes = minutes),
          approval_delay = 60 * 60 * 24,
          annotation = None,
          questions = None,
          qualifications = qualifications,
          layout_params = None,
          response_groups = None
          )
      return hit

'''
  Analysis Helpers
  -----

  General functions for checking experimental results
'''

class Analysis():
    def __init__(self):
      self.error_aggreements = {
          'true_positives' : [],
          'false_positives': [],
          'false_negatives': [] }


#
#     def flatten(self, iterables):
#       return (elem for iterable in iterables for elem in iterable)
#
#     def process_annotations(self, user=None, document=None, view=None, strict=False):
#       '''
#         This function returns all the dictionary annotations for a particular user and document or
#         for a given View Model
#       '''
#       if user is not None and document is not None:
#         annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user = user).all()
#       elif view is not None:
#         annotations = db.session.query(Annotation).filter_by(document = view.document).filter_by(user = view.user).all()
#       else:
#         raise ValueError("Not enough information given to retrieve annotations")
#
#       if strict and len(annotations) is 0:
#         raise ValueError( "No annotations available for this document")
#
#       # Convert all of the Anotation SQLAlchemy models to cleaned up dictionaries with
#       # cleaned whitespace and offsets
#       annotations = [ann.compare_view() for ann in annotations]
#       uniq_annotations = [dict(y) for y in set(tuple(x.items()) for x in annotations)]
#
#       # if len(annotations) is not len(uniq_annotations):
#       #   print annotations
#       #   print uniq_annotations
#       #   print ' - - - - - - - - - - '
#
#       return uniq_annotations
#

#
#     def determine_f(self, true_positive, false_positive, false_negative):
#       if true_positive + false_positive is 0:
#         return (0,0,0)
#
#       precision = true_positive / float(true_positive + false_positive)
#       recall = true_positive / float(true_positive + false_negative)
#
#       # print "Precision: {} Recall: {}".format(precision, recall)
#
#       if precision + recall > 0.0:
#         f = ( 2 * precision * recall ) / ( precision + recall )
#         return (precision, recall, f)
#       else:
#         return (0,0,0)
#
#     def figure_two(self, documents, match=0, experiment=2):
#       '''
#         Generates a table of NCBO vs K1-5 users annotations for the list of input Documents
#
#         Input:
#           documents - List of SQLAlchemy Document models
#           match - int for type of matching algorithmn to be used. 0 = Exact, 1 = Partial
#         Output: none (prints ASCII Table)
#       '''
#       # Setup the results dictionary witch will store all the annotations we're comparing
#       results = {'ncbo':[]}
#       for item in range(0,6): results[item] = []
#
#       for document in documents:
#         # Collect the list of Annotations models for the Golden Master and NCBO Annotator to use throughout
#         gm_annotations = self.process_annotations(user=User.query.get(2), document=document)
#         ncbo_annotations = self.process_annotations(user=User.query.get(1), document=document)
#
#         ncbo_score = self.calc_score(ncbo_annotations, gm_annotations, match)
#         results['ncbo'].append( ncbo_score )
#
#         # Collect all of the MTurk workers that viewed this document and assemble a list of all
#         # the annotations across the MTurk Workers for this document in the experiment
#         worker_views = self.get_experiment_hits(experiment, document)
#         workers_culmulative = self.get_workers_culmulative_annotations(worker_views)
#
#         # Build dictionary to store different K scores
#         k = {}
#         for item in range(0,8): k[item] = [] # This is higher b/c some anns might have higher agreements than # of users
#         # Looping through all the unique annotations to get their counts to actual
#         # submitted annotations for the workers results
#         uniq_anns = [dict(y) for y in set(tuple(x.items()) for x in workers_culmulative)]
#         for ann in uniq_anns:
#           # Put that annotation into the k for the # that it matches (how many times did workers agree on that
#           # particular annotation) and everything below it
#           # Ex: If an annotation matches 3 times, it also matches 2 times
#           for item in range(0, workers_culmulative.count(ann)): k[ item ].append( ann )
#
#         # Put the document K score annotations and append their TP/FP/FN counts to the K results
#         for i in range(0, 5): results[i].append( self.calc_score(k[i], gm_annotations, match) )
#
#       # We've now built up the results dictionary for our K scores and NCBO annotator for all the documents.
#       # Sum all the scores up, calculate their P/R/F and print it out
#       print "\t".join(['user ', 'p', 'r', 'f'])
#       for group in range(0,5):
#         results[group] = map(sum,zip(*results[group]))
#         results[group] = self.determine_f( results[group][0], results[group][1], results[group][2] )
#         print "\t".join(["{} ".format(group), "%.2f"%results[group][0], "%.2f"%results[group][1], "%.2f"%results[group][2]])
#
#     def get_workers_culmulative_annotations(self, worker_views):
#         # Collect all of the MTurk workers that viewed this document and assemble a list of all
#         # the annotations across the MTurk Workers for this document in the experiment
#         workers_culmulative = []
#         for worker_view in worker_views:
#           annotations = self.process_annotations(view=worker_view)
#           workers_culmulative.append( annotations )
#         # Flatten the list so we just have an array of all the combined workers annotations
#         workers_culmulative = [item for sublist in workers_culmulative for item in sublist]
#         return workers_culmulative
#


#
#     def util_annotation_length(self, experiment=2):
#       annotations = self.get_experiment_annotations(experiment)
#       annotations = [ann.compare_view() for ann in annotations]
#       annotations = [len(ann['text']) for ann in annotations]
#       total = float(len(annotations))
#       annotations = collections.Counter( annotations )
#
#       print "\t".join(["Annotation Length", "Occurances", "Percentage"])
#       for ann_len in annotations.items():
#         print "\t".join([str(ann_len[0]), str(ann_len[1]), "%.2f"%(ann_len[1]/total)])
#
#     def util_worker_contribution_counts(self, experiment=2):
#       hits = self.get_experiment_hits(experiment)
#
#       total = float(len(hits))
#
#       users = [hit.user.username for hit in hits]
#       users = collections.Counter( users )
#       users = dict((str(k), v) for k, v in users.iteritems())
#       print "\t".join(["User", "Submissions", "Percentage"])
#       for user in users:
#         print "\t".join([user, str(users[user]), "%.2f"%(users[user]/total) ])
#




#     def util_demographic(self, experiment=2):
#       annotations = self.get_experiment_annotations(experiment)
#       ips = [ann.player_ip for ann in annotations]
#       total = float(len(annotations))
#       ips = collections.Counter( ips )
#
#       print "\t".join(["IP", "Occurances", "Percentage", "City", "Country"])
#       for location in ips.items():
#         percent = (location[1]/ float(total)) * 100
#         r = requests.get('http://api.hostip.info/get_json.php?ip='+ location[0] +'&position=true').json()
#         print "\t".join([location[0], str(location[1]), "%.2f"%(location[1]/total), str(r.get('city', "None")), str(r.get('country_name', "None"))])
#
#       return True



#     def util_global_score(self, documents, experiment=2):
#       '''
#         Calculates the fp/fn/tp for a selection of user submissions
#       '''
#       for document in documents:
#         gm_annotations = self.process_annotations(user = User.query.get(2), document = document)
#
#         # Collect all of the MTurk workers that viewed this document and assemble a list of all
#         # the annotations across the MTurk Workers for this document in the experiment
#         worker_views = self.get_experiment_hits(experiment, document)
#         workers_culmulative = self.get_workers_culmulative_annotations( worker_views )
#         uniq_workers_culmulative = [dict(y) for y in set(tuple(x.items()) for x in workers_culmulative)]
#
#         # Runs the comparision between the workers and the gold master, saving to our class vars in the process
#         self.calc_score(uniq_workers_culmulative, gm_annotations)
#
#
#       shared_keys = []
#       for key in self.error_aggreements.keys():
#         self.error_aggreements[key] = collections.Counter(self.error_aggreements[key])
#         self.error_aggreements[key] = dict((str(k), v) for k, v in self.error_aggreements[key].iteritems())
#         shared_keys.append( self.error_aggreements[key].keys() )
#
#       shared_keys = list(set( self.flatten(shared_keys) ))
#       shared_keys.sort()
#       for key in shared_keys:
#         tp = str( self.error_aggreements['true_positives'].get(key, 0) )
#         fp = str( self.error_aggreements['false_positives'].get(key, 0) )
#         fn = str( self.error_aggreements['false_negatives'].get(key, 0) )
#         print "\t".join([key, tp, fp, fn])
#
