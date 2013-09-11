from flask import jsonify
from flask.ext.restful import reqparse, Resource

from models import db, User, Annotation

import re

class Gold(Resource):
  def get(self, doc_id):
    # Get the gold bot account to make the db entries
    user = db.session.query(User).get(2)
    annotations = db.session.query(Annotation).filter_by(document_id = doc_id).filter_by(user = user).all()

    return jsonify(objects=[ann.json_view() for ann in annotations])
