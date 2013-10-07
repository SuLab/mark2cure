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

document_parser = reqparse.RequestParser()
document_parser.add_argument('document_id',   type=int,   location='json')
document_parser.add_argument('title',         type=str,   location='json')
document_parser.add_argument('annotations',   type=list,  location='json')
# document_parser.add_argument('q',             type=str)

# Will only ever GET from Documents to initialize
# Will only ever PUT to Documents to add annotations
# Will serve 'heatmap' models in fetch based on word index sums

class Documents(Resource):
    def get(self, doc_id=None):
        print "docss"
        print current_user

        if doc_id:
          document = db.session.query(Document).get(doc_id)
          return jsonify(objects=[document.json_view(current_user)])
        else:
          documents = Document.query.limit(30).all()
          return jsonify(objects=[i.json_view(current_user) for i in documents])

    def put(self, doc_id):
        args = document_parser.parse_args()
        user = current_user
        document = Document.query.get(doc_id)
        env =  request.environ

        view = View(user, document);
        db.session.add(view)
        db.session.commit()

        for ann in args['annotations']:
          ann = Annotation( ann['kind'],
                            ann['type'],
                            ann['text'],
                            ann['start'],
                            ann['length'],
                            ann['stop'],
                            user,
                            document,
                            env.get('HTTP_USER_AGENT'),
                            env.get('REMOTE_ADDR'),
                            None
                          );
          db.session.add(ann)

        db.session.commit()
        return args, 201


