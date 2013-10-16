from flask import request, jsonify
from flask.ext.restful import reqparse, Resource
from flask_login import login_user, current_user

# from sqlalchemy import *
# from sqlalchemy.orm import *
# from sqlalchemy.ext.declarative import declarative_base
# from sqlalchemy.sql.expression import ClauseElement
# from collections import OrderedDict

from ..models import User, Document, Annotation, View
from ..core import db
from mark2cure.settings import *

document_parser = reqparse.RequestParser()
document_parser.add_argument('document_id',   type=int,   location='json')
document_parser.add_argument('title',         type=str,   location='json')
document_parser.add_argument('annotations',   type=list,  location='json')

# Will only ever GET from Documents to initialize
# Will only ever PUT to Documents to add annotations
# Will serve 'heatmap' models in fetch based on word index sums

class Documents(Resource):
    def get(self, doc_id=None):
        if doc_id:
          document = db.session.query(Document).get(doc_id)
          return jsonify(objects=[document.json_view(current_user)])
        else:
          documents = Document.query.limit(30).all()
          return jsonify(objects=[i.json_view(current_user) for i in documents])

    def put(self, doc_id):
        args = document_parser.parse_args()
        document = Document.query.get(doc_id)
        env =  request.environ

        view = View(current_user, document);
        db.session.add(view)
        db.session.commit()

        for ann in args['annotations']:
          ann = Annotation( ann['kind'],
                            ann['type'],
                            ann['text'],
                            ann['start'],
                            ann['length'],
                            ann['stop'],
                            current_user,
                            document,
                            env.get('HTTP_USER_AGENT'),
                            env.get('REMOTE_ADDR'),
                            None
                          );
          db.session.add(ann)
        db.session.commit()

        # Check document and mturk status
        if document.validate and current_user.mturk:
          # This is a document that requires validation
          user_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user = current_user).all()
          gold_annotations = db.session.query(Annotation).filter_by(document = document).filter_by(user_id = 2).all()

          user_annotations = [ann.compare_view() for ann in user_annotations]
          gold_annotations = [ann.compare_view() for ann in gold_annotations]

          user_matches = len([ann for ann in user_annotations if ann in gold_annotations])

          print user_matches, ( len(user_annotations) - user_matches), len(user_annotations), len(gold_annotations)

          # mtc = MTurkConnection( aws_access_key_id = AWS_ACCESS_ID,
          #                        aws_secret_access_key = AWS_SECRET_KEY,
          #                        host = AWS_HOST)

          # score = int(mtc.get_qualification_score(AWS_QUAL_GM_SCORE, worker))
          # if(len(user_matches) > 0):
          #   score += 1
          # else:
          #   score -= 1

          # mtc.update_qualification_score(AWS_QUAL_GM_SCORE, worker, score)

        return args, 201

