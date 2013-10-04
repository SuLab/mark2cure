from ..core import db
import datetime

class Message(db.Model):
    id        = db.Column(db.Integer, primary_key=True)
    message   = db.Column(db.Text)
    created   = db.Column(db.DateTime)

    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))

    def __init__(self, message, user):
        self.message    = message
        self.created    = datetime.datetime.utcnow()
        self.user       = user
