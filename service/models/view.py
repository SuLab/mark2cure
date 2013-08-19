from . import db
import datetime

class View(db.Model):
    id          = db.Column(db.Integer, primary_key=True)
    created     = db.Column(db.DateTime)

    user_id       = db.Column(db.Integer, db.ForeignKey('user.id'))
    document_id   = db.Column(db.Integer, db.ForeignKey('document.id'))

    def __init__(self, user, document):
        self.created    = datetime.datetime.utcnow()
        self.user       = user
        self.document   = document
