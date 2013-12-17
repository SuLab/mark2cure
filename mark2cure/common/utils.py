from django.conf import settings
from datetime import datetime
from django.contrib.auth.models import User

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
      return ExternalQuestion("https://mark2cure.org/document/"+str(doc_id), 800)

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
                                        '<h2>Instructions</h2>'
                                        '<ol>'
                                          '<li>'
                                            '<h3>Highlight <u>all</u> diseases and disease abbreviations</h3>'
                                            '<img alt="Highlight all diseases and disease abbreviations" src="http://mark2cure.org/static/images/experiment/3/1.png" />'
                                            '<img alt="Highlight all diseases and disease abbreviations" src="http://mark2cure.org/static/images/experiment/3/2.png" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                            '<h3>Highlight the longest span of text specific to a disease</h3>'
                                            '<img alt="Highlight the longest span of text specific to a disease" src="http://mark2cure.org/static/images/experiment/3/3.png" />'
                                            '<img alt="Highlight the longest span of text specific to a disease" src="http://mark2cure.org/static/images/experiment/3/4.png" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                            '<h3>Highlight disease conjunctions as single, long spans.</h3>'
                                            '<img alt="Highlight disease conjunctions as single, long spans" src="http://mark2cure.org/static/images/experiment/3/5.png" />'
                                            '<img alt="Highlight disease conjunctions as single, long spans" src="http://mark2cure.org/static/images/experiment/3/6.png" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                            '<h3>Highlight symptoms - physical results of having a disease</h3>'
                                            '<img alt="Highlight symptoms - physical results of having a disease" src="http://mark2cure.org/static/images/experiment/3/7.png" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                            '<h3>Highlight <u>all</u> occurrences of disease terms</h3>'
                                            '<img alt="Highlight all occurrences of disease terms" src="http://mark2cure.org/static/images/experiment/3/8.png" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                          '<li>'
                                            '<h3>Highlight <u>all</u> diseases, disease groups and key disease symptoms</h3>'
                                            '<img alt="Highlight all diseases, disease groups and key disease symptoms" src="http://mark2cure.org/static/images/experiment/3/9.png" />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                            '<br />'
                                          '</li>'
                                        '</ol>'))

      #
      # Questions to ask the Worker
      #
      instructions = "Select all and only the terms that should be highlighted for each text segment (don't select terms that overlap with each other in the text):"
      paragraph1 = "Myotonic dystrophy ( DM ) is associated with a ( CTG ) n trinucleotide repeat expansion in the 3-untranslated region of a protein kinase-encoding gene , DMPK , which maps to chromosome 19q13 . 3 . "
      paragraph2 = "Germline mutations in BRCA1 are responsible for most cases of inherited breast and ovarian cancer . However , the function of the BRCA1 protein has remained elusive . As a regulated secretory protein , BRCA1 appears to function by a mechanism not previously described for tumour suppressor gene products."
      paragraph3 = "We report about Dr . Kniest , who first described the condition in 1952 , and his patient , who , at the age of 50 years is severely handicapped with short stature , restricted joint mobility , and blindness but is mentally alert and leads an active life .  This is in accordance with molecular findings in other patients with Kniest dysplasia and..."

      # # # # # # # # #
      #
      # Paragraph 1
      #
      # # # # # # # # #
      #abcdefghijklmnopqrstuvwxyz
      qc = QuestionContent()
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
      qc.append_field('Text', 'DMPK')
      q7 = Question(identifier = 'term_selection_7', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t7"), ("False", "f7")])))

      # # # # # # # # #
      #
      # Paragraph 2
      #
      # # # # # # # # #
      qc = QuestionContent()
      qc.append_field('Title', paragraph2)
      qc.append_field('Text', 'Germline mutations')
      q8 = Question(identifier = 'term_selection_8', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t8"), ("False", "f8")])))
      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'inherited breast and ovarian cancer')
      q9 = Question(identifier = 'term_selection_9', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t9"), ("False", "f9")])))
      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'breast')
      q10 = Question(identifier = 'term_selection_10', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t10"), ("False", "f10")])))
      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'ovarian cancer')
      q11 = Question(identifier = 'term_selection_11', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t11"), ("False", "f11")])))
      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'cancer')
      q12 = Question(identifier = 'term_selection_12', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t12"), ("False", "f12")])))
      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'tumour')
      q13 = Question(identifier = 'term_selection_13', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t13"), ("False", "f13")])))
      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'tumour suppressor')
      q14 = Question(identifier = 'term_selection_14', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t14"), ("False", "f14")])))


      # # # # # # # # #
      #
      # Paragraph 3
      #
      # # # # # # # # #
      qc = QuestionContent()
      qc.append_field('Title', paragraph3)
      qc.append_field('Text', 'age of 50 years')
      q15 = Question(identifier = 'term_selection_15', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t15"), ("False", "f15")])))
      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'short stature')
      q16 = Question(identifier = 'term_selection_16', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t16"), ("False", "f16")])))
      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'restricted joint mobility')
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
      qc.append_field('Text', 'blindness')
      q20 = Question(identifier = 'term_selection_20', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t20"), ("False", "f20")])))
      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'dysplasia')
      q21 = Question(identifier = 'term_selection_21', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t21"), ("False", "f21")])))
      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'Kniest dysplasia')
      q22 = Question(identifier = 'term_selection_22', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t22"), ("False", "f22")])))
      qc = QuestionContent()
      qc.append_field('Title', "")
      qc.append_field('Text', 'molecular findings')
      q23 = Question(identifier = 'term_selection_23', content = qc, is_required = True,
          answer_spec = AnswerSpecification(SelectionAnswer(
            style='radiobutton',
            selections=[("True", "t23"), ("False", "f23")])))


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
                              <SelectionIdentifier>t9</SelectionIdentifier>
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
                              <SelectionIdentifier>f15</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_16</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>t16</SelectionIdentifier>
                              <AnswerScore>1</AnswerScore>
                            </AnswerOption>
                          </Question>

                          <Question>
                          <QuestionIdentifier>term_selection_17</QuestionIdentifier>
                            <AnswerOption>
                              <SelectionIdentifier>t17</SelectionIdentifier>
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
                              <SelectionIdentifier>f21</SelectionIdentifier>
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

                        </AnswerKey>'''

      qual_test = self.mtc.update_qualification_type(settings.AWS_QUAL_TEST_3,
      # qual_test = self.mtc.create_qualification_type(
        # name = 'TEST3 :: Annotation Instructions & Qualification Questions',
        description = 'Detailed annotation instructions. Multiple-choice questions to assess concept understanding.',
        status = 'Active',
        test = question_form,
        answer_key = answer_logic,
        retry_delay = 1,
        test_duration = 20 * 60)

      return qual_test

  # Actionable methods
  def hit_for_document(self, doc_id, max_assignments = 5, reward = 0.06, minutes = 4, title="Highlight diseases in paragraph"):
      description = ('You will be presented with text from the biomedical literature which we believe may help resolve some important medical related questions. The task is to highlight words and phrases in that text which are, or are highly related to, diseases.  This work will help advance research in cancer and many other diseases!')
      keywords = 'science, abstract, primary literature, annotation, disease, text, highlight, annotation, medicine, term recognition'

      qualifications = Qualifications()
      # Add the step instructions and basic test
      qualifications.add( Requirement(settings.AWS_QUAL_TEST_3, "GreaterThanOrEqualTo", 20) )

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
