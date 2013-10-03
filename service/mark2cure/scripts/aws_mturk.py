from boto.mturk.connection import MTurkConnection
from boto.mturk.question import QuestionContent, ExternalQuestion


mtc = MTurkConnection(aws_access_key_id=ACCESS_ID,
                      aws_secret_access_key=SECRET_KEY,
                      host=HOST)

title = 'Help annotate scientific articles'
description = ('Visit a website and highlight diseases that are present in a paragraph.')
keywords = 'science, annotation, disease'

#---------------  BUILD OVERVIEW -------------------
 
overview = Overview()
overview.append_field('Title', 'Give your opinion on this website')
overview.append(FormattedContent('<a target="_blank"'
                                 ' href="http://www.toforge.com">'
                                 ' Mauro Rocco Personal Forge</a>'))
 
#---------------  BUILD QUESTION 1 -------------------
 
qc1 = QuestionContent()
qc1.append_field('Title','How looks the design ?')
 
fta1 = SelectionAnswer(min=1, max=1,style='dropdown',
                      selections=ratings,
                      type='text',
                      other=False)
 
q1 = Question(identifier='design',
              content=qc1,
              answer_spec=AnswerSpecification(fta1),
              is_required=True)
 
#---------------  BUILD QUESTION 2 -------------------
 
qc2 = QuestionContent()
qc2.append_field('Title','Your personal comments')
 
fta2 = FreeTextAnswer()
 
q2 = Question(identifier="comments",
              content=qc2,
              answer_spec=AnswerSpecification(fta2))
 
#--------------- BUILD THE QUESTION FORM -------------------
 
question_form = QuestionForm()
question_form.append(overview)
question_form.append(q1)
question_form.append(q2)
 
#--------------- CREATE THE HIT -------------------
 
mtc.create_hit(questions=question_form,
               max_assignments=1,
               title=title,
               description=description,
               keywords=keywords,
               duration = 60*5,
               reward=0.05)
