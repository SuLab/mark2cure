from flask import jsonify
from flask.ext.restful import reqparse, Resource
from flask_login import current_user

from ..core import db
from ..models import User, Annotation

import re

def clean(text):
  word = re.sub(r'\W+', ' ', text).lower().strip()
  return word[:-1] if word in ['cancers', 'chordomas'] else word

class Network(Resource):
  def get(self):
    #
    # THIS ALGO IMPLEMENTS ALL USER ANNOTATIONS ACCROSS ALL DOCS
    #
    node_list = []
    link_list = []

    # (TODO) Refine to X # of articles, X days ago, ...
    annotations = db.session.query(Annotation).filter_by(user = current_user).all()

    # Saving the unique arrays of annotations and documents
    ann_arr = list(set( [clean(ann.text)  for ann in annotations] ))
    doc_arr = list(set( [ann.document     for ann in annotations] ))

    # Put the uniq docs at the top of the list then uniq annotation terms after the documents
    [node_list.append({"name": doc.title, "group": 1})  for doc in doc_arr]
    [node_list.append({"name": ann, "group": 2})        for ann in ann_arr]

    # For each of the documents with annotations
    for doc_idx, doc in enumerate( doc_arr ):
      # Collect all of the cleaned annotations for that document by the current user
      doc_ann = [clean(ann.text) for ann in doc.annotations if ann.user.id is current_user.id]

      for ann in list(set(doc_ann)):
        # For this document, connect the doc to the annotation and weight appropriately
        link_list.append({"source"  : doc_idx,
                          "target"  : len(doc_arr) + ann_arr.index(ann),
                          "value"   : doc_ann.count(ann) })

    return jsonify(nodes=node_list, links=link_list)
