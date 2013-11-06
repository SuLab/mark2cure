from flask import request, jsonify
from flask.ext.restful import reqparse, Resource
from flask_login import login_user, current_user
from flask_mail import Message

# from sqlalchemy import *
# from sqlalchemy.orm import *
# from sqlalchemy.ext.declarative import declarative_base
# from sqlalchemy.sql.expression import ClauseElement
# from collections import OrderedDict
from sqlalchemy import desc

from ..models import User, Document, Annotation, View
from ..core import db, mail
from mark2cure.settings import *
from mark2cure.manage.analysis import gold_matches
from mark2cure.manage.aws import Turk

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

        views = db.session.query(View).filter_by(user = current_user).filter_by( document = document ).all()
        if len(views):
          if current_user.mturk:
            t = Turk()
            t.mtc.block_worker(current_user.username, "Attempted to submit same document multiple times.")
          raise ValueError("Cannot submit a document twice")
        else:
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

          # This makes it easier to track in the DB later on
          if(current_user.mturk):
            ann.experiment = 2

          db.session.add(ann)
        db.session.commit()

        if document.validate and current_user.mturk:
            # If this is a validate document, check the user's history, if it's their 3rd submission
            # or more run test to potentially fail if poor performance
            valid_views = db.session.query(View).filter_by(user = current_user).filter( View.document.has(validate=1) ).order_by( desc(View.created) ).limit(3).all()
            if len(valid_views) is 3:
              if sum(1 for x in valid_views if gold_matches(x.user, x.document) >= 1) is not 3:
                print "failed"
                t = Turk()
                t.mtc.block_worker(current_user.username, "Failed to properly answer golden master performance documents")

        return args, 201

