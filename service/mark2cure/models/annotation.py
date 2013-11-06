from ..core import db
import datetime

class Annotation(db.Model):
    id          = db.Column(db.Integer,   primary_key=True)

    kind    = db.Column(db.Integer)
    type    = db.Column(db.String(12))
    text    = db.Column(db.Text)
    start   = db.Column(db.Integer)
    length  = db.Column(db.Integer)
    stop    = db.Column(db.Integer)

    created = db.Column(db.DateTime)

    user_id       = db.Column(db.Integer, db.ForeignKey('user.id'))
    document_id   = db.Column(db.Integer, db.ForeignKey('document.id'))
    concept_id    = db.Column(db.Integer, db.ForeignKey('concept.id'))

    user_agent  = db.Column(db.String(150))
    player_ip   = db.Column(db.String(30))
    experiment = db.Column(db.Integer)

    def __init__(self, kind, type, text, start, length, stop, user, document, ua, pip, concept):
        self.kind     = kind
        self.type     = type
        self.text     = text
        self.start    = start
        self.length   = length
        self.stop     = stop

        self.created  = datetime.datetime.utcnow()

        self.user     = user
        self.document = document
        self.concept  = concept

        self.user_agent = ua
        self.player_ip = pip

    # Required for administrative interface
    def __unicode__(self):
        return self.text

    def __repr__(self):
        return '<Ann %r>' % self.type

    def json_view(self):
      return {  'id'        : self.id,
                'kind'      : self.kind,
                'type'      : self.type,
                'text'      : self.text,
                'start'     : self.start,
                'length'    : self.length,
                'stop'      : self.stop,
                'created'   : self.created.isoformat() }

    def compare_view(self):
      # Returns back the text dictionary for comparision
      offset = len(self.text) - len(self.text.lstrip())
      return {
                'text'      : self.text.strip().lower(),
                'start'     : int(self.start)+offset,
              }
