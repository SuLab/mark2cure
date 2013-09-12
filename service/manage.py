from flask import Flask
from flask.ext.script import Manager, Command, Option

from models import *

import settings, requests, re, csv, datetime
import xml.etree.ElementTree as ET
# from Bio import Entrez, Medline
from bs4 import BeautifulSoup, NavigableString

app = Flask(__name__)
# configure your app

app.config['SQLALCHEMY_DATABASE_URI'] = settings.DB_URI
db.init_app(app)

manager = Manager(app)

class Create(Command):
    "Create a database if doesn't already exist"

    def run(self):
      db.create_all()
      print ":)"

class Heatmap(Command):
    "Calculates the heatmap array"

    def run(self):
      documents = Document.query.all()
      for document in documents:
        print "Running document {0}".format( document.document_id )
        # Select all the annotations for this document
        annotations = db.session.query(Annotation).filter_by(document = document)
        step = 0
        length = 0
        pop_arr = []

        # Iterate over each word in the text
        for idx, word in enumerate(document.text.split()):
          length = len(word)
          step = step + length + 1;
          wordstart = step - length - 1,

          # Count which words have annotations in the db for that doc, across all users
          # Number of annotations that encapsulate that word
          count = 0
          for ann in annotations:
            if wordstart[0] >= ann.start and wordstart[0] <= (ann.start+ann.length):
              count = count + 1
          pop_arr.append(count);

        document.cache = ', '.join([str(x) for x in pop_arr])

        # Save every document instead of once incase some doc crashes
        db.session.commit()
      print "Complete"

class Annotate(Command):
    "Populate the documents with the NCBI Annotator"

    def run(self):
      # Get the annotator bot account to make the db entries
      user = db.session.query(User).get(1)

      # NCBO Annotator http request
      payload = { 'apikey'                    : settings.NCBO_API_KEY,
                  'stopWords'                 : settings.STOP_WORDS,
                  'minTermSize'               : 5,
                  'withSynonyms'              : 'true',
                  # 1009 = Human disease ontology
                  'ontologiesToKeepInResult'  : '1009',
                  'isVirtualOntologyId'       : 'true',
                  # T047 -- Disease or Syndrome
                  # 'semanticTypes'             : 'T047',
                  'textToAnnotate'            : '',
                }

      # Go over all the documents
      documents = Document.query.all()
      for document in documents:
        # If the current document doesn't have at least 4 annotations from the bot, try to get more...
        if db.session.query(Annotation).filter_by(document = document).filter_by(user = user).count() <= 3:
          print "Running document {0}".format( document.document_id )
          payload['textToAnnotate'] = document.text
          r = requests.post("http://rest.bioontology.org/obs/annotator", data=payload)

          if r.ok:
            try:
              root = ET.fromstring( r.text.decode('utf-8') )
              # Make array of all relevant NCBO results
              for ann in root.iter('annotationBean'):
                ctx = ann.find('context')
                start = int(ctx.find('from').text)-1
                stop = int(ctx.find('to').text)
                annotation = document.text[start:stop]
                concept_url =  ctx.find('term').find('concept').find('fullId').text

                concept = db.session.query(Concept).filter_by(concept_id = concept_url).first()
                if concept is None:
                  concept = Concept(concept_url)
                  db.session.add(concept)
                  db.session.commit()

                # Find those results in the original abstract we had
                ann = Annotation( 0,
                                  'disease',
                                  annotation,
                                  start,
                                  len(annotation),
                                  stop,
                                  user,
                                  document,
                                  'annotator_bot_1.0',
                                  '',
                                  concept
                                );
                db.session.add(ann)

              # Save every document instead of once incase some doc crashes
              db.session.commit()
            except Exception:
              pass
      print ":)"

# class Import(Command):
#     MAX_COUNT = 100
#     Entrez.email = 'nanis@scripps.edu'
#
#     for term in ['chordoma', 'cancer']:
#       h = Entrez.esearch(db='pubmed', retmax=MAX_COUNT, term=term)
#       result = Entrez.read(h)
#       ids = result['IdList']
#       h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
#       records = Medline.parse(h)
#
#       #
#       # Reference to abbreviations: http://www.nlm.nih.gov/bsd/mms/medlineelements.html
#       #
#       for record in records:
#         if record.get('TI') and record.get('AB') and record.get('PMID') and record.get('CRDT'):
#           instance = db.session.query(Document).filter_by(document_id = record.get('PMID')).count()
#           if instance == 0:
#             doc = Document( record.get('PMID'),
#                             record.get('AB'),
#                             record.get('TI'),
#                             datetime.datetime.strptime(record.get('CRDT')[0], '%Y/%m/%d %H:%M')
#                           )
#             db.session.add(doc)
#           db.session.commit()
    # print "Import"

def strip_tags(html, invalid_tags):
    soup = BeautifulSoup(html)

    for tag in soup.findAll(True):
        if tag.name in invalid_tags:
            s = ""

            for c in tag.contents:
                if not isinstance(c, NavigableString):
                    c = strip_tags(unicode(c), invalid_tags)
                s += unicode(c)

            tag.replaceWith(s)

    return soup

class SolidGold(Command):
    "Populate the documents with the NCBI SolidGold"

    option_list = (
        Option('--path', '-p', dest='path'),
    )

    def run(self, path):
      # Get the gold bot account to make the db entries
      user = db.session.query(User).get(2)
      with open('../assets/NCBI_corpus/'+ path +'.txt','r') as f:
        reader = csv.reader(f, delimiter='\t')
        for num, title, text in reader:
          text = text.replace('<category=', '<category type=')
          plain_text = str(strip_tags(text, ['category']))

          title = title.replace('<category=', '<category type=')
          plain_title = str(strip_tags(title, ['category']))

          soup = BeautifulSoup(text)

          # Add the document to the database
          doc = Document( int(num),
                          plain_text,
                          plain_title,
                          datetime.datetime.utcnow(),
                          path )
          db.session.add(doc)
          db.session.commit()

          cats = set([cat.text for cat in soup.findAll('category')])
          for cat in cats:
            for m in re.finditer( cat, plain_text ):
              # print( cat.text, cat['type'], m.start(), m.end() )
              ann = Annotation( 0,
                                'disease',
                                cat,
                                m.start(),
                                len(cat),
                                m.end(),
                                user,
                                doc,
                                'gold_1.0',
                                '',
                                None );
              db.session.add(ann)
              # Save every document instead of once incase some doc crashes
              db.session.commit()

class Compare(Command):
    "F Score"

    def run(self):
      document = db.session.query(Document).get(28)
      annotator_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user_id = 1).all()
      gold_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user_id = 2).all()

      print [ann.compare_view() for ann in annotator_annotations]
      print
      print [ann.compare_view() for ann in gold_annotations]

manager.add_command('heatmap', Heatmap())
manager.add_command('annotate', Annotate())
manager.add_command('create', Create())
manager.add_command('gold', SolidGold())
manager.add_command('compare', Compare())

if __name__ == "__main__":
    manager.run()
