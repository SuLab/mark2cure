from ..core import db

from .view import View
from .annotation import Annotation

import datetime

class Document(db.Model):
    id          = db.Column(db.Integer,   primary_key=True)
    document_id = db.Column(db.Integer)
    text        = db.Column(db.Text)
    title       = db.Column(db.Text)

    created     = db.Column(db.DateTime)
    cache       = db.Column(db.Text)
    source      = db.Column(db.String(200))

    views       = db.relationship('View',             backref=db.backref('document',  lazy='select'))
    annotations = db.relationship('Annotation',       backref=db.backref('document',  lazy='select'))
    quest_relations = db.relationship('QuestRelation', backref=db.backref('document',  lazy='select'))

    def __init__(self, document_id, text, title, created, source):
        self.document_id  = document_id
        self.text         = text
        self.title        = title
        if created is None:
          self.created      = datetime.datetime.utcnow()
        self.created = created
        self.source = source

    # Required for administrative interface
    def __unicode__(self):
        return self.title

    def __repr__(self):
        return '<Document %r>' % self.document_id

    def json_view(self, user):
        print user

        if user.is_anonymous():
          viewed = False
          annotations = []
        else:
          viewed = True if db.session.query(View).filter_by(user = user).filter_by(document = self).first() else False
          annotations = db.session.query(Annotation).filter_by(user = user).filter_by(document = self).all()

        popularity = self.cache.split(", ") if self.cache else ["0"]

        return {  'id'          : self.id,
                  'document_id' : self.document_id,
                  'text'        : self.text,
                  'title'       : self.title,
                  'created'     : self.created.isoformat(),
                  # relationships
                  'complete'    : viewed,
                  'annotations' : [i.json_view() for i in annotations],
                  'popularity'  : popularity
                  }
