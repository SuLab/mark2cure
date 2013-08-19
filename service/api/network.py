from flask import jsonify
from flask.ext.restful import reqparse, Resource

from models import db, User, Annotation

import re

def clean(text):
  word = re.sub(r'\W+', '', text).lower().strip()
  return word[:-1] if word in ['cancers', 'chordomas'] else word

network_parser = reqparse.RequestParser()
network_parser.add_argument('api_key',       type=str,   location='cookies')
class Network(Resource):
  def get(self):
    args = network_parser.parse_args();
    user = db.session.query(User).filter_by(api_key = args['api_key']).first()

    node_list = []
    link_list = []
    #
    # THIS ALGO IMPLEMENTS ALL SHARED ANNOTATIONS ACCROSS ALL DOCS
    #
    annotations = db.session.query(Annotation).filter_by(user = user).all()
    ann_arr = []
    doc_arr = []
    for ann in annotations:
      ann_text = clean(ann.text)
      ann_arr.append( ann_text )
      doc_arr.append( ann.document )

    # Saving the clean sets
    ann_arr = list(set(ann_arr))
    doc_arr = list(set(doc_arr))

    # Put the uniq docs at the top of the list
    for doc in doc_arr:
      node_list.append({"name": doc.title, "group": 1})

    # Put the uniq annotation terms after the documents
    for ann in ann_arr:
      node_list.append({"name": ann, "group": 2})

    # For each of the documents with annotations
    for doc_idx, doc in enumerate( doc_arr ):
      # Collect all of the cleaned annotations for that document
      doc_ann = [clean(a.text) for a in doc.annotations] 
      doc_ann_uniq = list(set(doc_ann))
      for ann in doc_ann_uniq:
        if(ann in ann_arr):
          # For this document, connect the doc to the annotation
          link_list.append({"source": doc_idx, "target": len(doc_arr) + ann_arr.index(ann) , "value": doc_ann.count(ann) })

    return jsonify(nodes=node_list, links=link_list)
