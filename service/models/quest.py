from . import db
import datetime

class Quest(db.Model):
    id          = db.Column(db.Integer, primary_key=True)
    name      = db.Column(db.String(200))
    created     = db.Column(db.DateTime)

    user_id       = db.Column(db.Integer, db.ForeignKey('user.id'))

    quest_relations = db.relationship('QuestRelation', backref=db.backref('quest',  lazy='select'))

    def __init__(self, name, user):
        self.name = name
        self.created    = datetime.datetime.utcnow()
        self.user       = user
