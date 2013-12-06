'''
  Calculates the heatmap array
'''

# def run(self):
#   for document in Document.query.all():
#     print "Running document {0}".format( document.document_id )
#
#
#     # Select all the annotations for this document
#     annotations = db.session.query(Annotation).filter_by(document = document)
#
#
#     step = 0
#     length = 0
#     pop_arr = []
#
#     # Iterate over each word in the text
#     for idx, word in enumerate(document.text.split()):
#       length = len(word)
#       step = step + length + 1;
#       wordstart = step - length - 1,
#
#       # Count which words have annotations in the db for that doc, across all users
#       # Number of annotations that encapsulate that word
#       count = 0
#       for ann in annotations:
#         if wordstart[0] >= ann.start and wordstart[0] <= (ann.start+ann.length):
#           count = count + 1
#       pop_arr.append(count);
#
#     document.cache = ', '.join([str(x) for x in pop_arr])
#
#     # Save every document instead of once incase some doc crashes
#     db.session.commit()
#   print "Complete"
#
#

'''
  Adds the NCBO annotations
'''

# import requests, re, csv, datetime
# import xml.etree.ElementTree as ET
#
# class Annotate(Command):
#     "Populate the documents with the NCBI Annotator"
#
#     def __init__(self):
#       # NCBO Annotator http request
#       self.payload = { 'apikey'                    : NCBO_API_KEY,
#                 'stopWords'                 : STOP_WORDS,
#                 'minTermSize'               : 5,
#                 'withSynonyms'              : 'true',
#                 # 1009 = Human disease ontology
#                 'ontologiesToKeepInResult'  : '1009',
#                 'isVirtualOntologyId'       : 'true',
#                 # T047 -- Disease or Syndrome
#                 # 'semanticTypes'             : 'T047',
#                 'textToAnnotate'            : '',
#               }
#
#     def handle_request(self, request, annotator, document):
#       if request.ok:
#         try:
#           root = ET.fromstring( request.text.decode('utf-8') )
#
#           # print request.text.decode('utf-8')
#           # Make array of all relevant NCBO results
#           for ann in root.iter('annotationBean'):
#             ctx = ann.find('context')
#             score = int(ann.find('score').text)
#             start = int(ctx.find('from').text)-1
#             stop = int(ctx.find('to').text)
#
#             annotation = document.text[start:stop]
#             concept_url =  ctx.find('term').find('concept').find('fullId').text
#
#             if score >= annotator.score:
#               concept = db.session.query(Concept).filter_by(concept_id = concept_url).first()
#               if concept is None:
#                 concept = Concept(concept_url)
#                 db.session.add(concept)
#                 db.session.commit()
#
#               # Find those results in the original abstract we had
#               ann = Annotation( 0,
#                                 'disease',
#                                 annotation,
#                                 start,
#                                 len(annotation),
#                                 stop,
#                                 annotator.user,
#                                 document,
#                                 annotator.user.username,
#                                 '',
#                                 concept
#                               );
#               db.session.add(ann)
#
#             # Save every document instead of once incase some doc crashes
#             db.session.commit()
#         except Exception:
#           pass
#
#     def update_annotations(self):
#       # Get the annotator bot account to make the db entries
#       annotators = db.session.query(Ncbo).all()
#
#       # Go over all the documents
#       # documents = Document.query.all()
#       documents = db.session.query(Document).\
#           filter_by(source = 'NCBI_corpus_development').\
#           all()
#
#       for document in documents:
#         for annotator in annotators:
#           # If the current document doesn't have at least 4 annotations from the bot, try to get more...
#           current_anns = db.session.query(Annotation).\
#               filter_by(document = document).\
#               filter_by( user = annotator.user ).\
#               count()
#
#           if current_anns <= 3:
#             print annotator.user.username, document.document_id
#
#             # Config settings to change per request
#             self.payload['textToAnnotate'] = document.text
#             self.payload['minTermSize'] = annotator.min_term_size
#
#             r = requests.post("http://rest.bioontology.org/obs/annotator", data=self.payload)
#             self.handle_request(r, annotator, document)
#             print '- - - - - - - - - - - - -'
#
#     def util_add_users(self):
#       for minsize in range(2,8):
#         for score in range(5, 100, 10):
#           ncbo = Ncbo(minsize, score)
#           db.session.add(ncbo)
#           db.session.commit()
#       return True
#
