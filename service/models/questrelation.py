from . import db
import datetime

class QuestRelation(db.Model):
    id          = db.Column(db.Integer, primary_key=True)
    created     = db.Column(db.DateTime)

    quest_id       = db.Column(db.Integer, db.ForeignKey('quest.id'))
    document_id   = db.Column(db.Integer, db.ForeignKey('document.id'))

    def __init__(self, quest, document):
        self.created    = datetime.datetime.utcnow()
        self.quest      = quest
        self.document   = document

    def json_view(self):
      return { 'document_id'  : self.document.id,
                'quest_id'    : self.quest.id,
                'created'     : self.created.isoformat()  }
