from . import db
import datetime

class User(db.Model):
    id          = db.Column(db.Integer,     primary_key=True)

    username    = db.Column(db.String(120), unique=True)
    email       = db.Column(db.String(120))
    experience  = db.Column(db.Integer())

    created     = db.Column(db.DateTime)

    feedback_0  = db.Column(db.Integer())
    feedback_1  = db.Column(db.Integer())
    feedback_2  = db.Column(db.Integer())
    feedback_3  = db.Column(db.Integer())

    first_run   = db.Column(db.Boolean())
    api_key     = db.Column(db.Text)

    views       = db.relationship('View',         backref=db.backref('user',  lazy='select'))
    messages    = db.relationship('Message',      backref=db.backref('user',  lazy='select'))
    annotations = db.relationship('Annotation',   backref=db.backref('user',  lazy='select'))

    # Required for administrative interface
    def __unicode__(self):
        return self.email

    def __init__(self, username, email, experience, api):
        self.username   = username
        self.email      = email
        self.experience = experience
        self.created    = datetime.datetime.utcnow()

        self.feedback_0 = -1
        self.feedback_1 = -1
        self.feedback_2 = -1
        self.feedback_3 = -1

        self.first_run  = True
        self.api_key    = api

    def __repr__(self):
        return '<User %r>' % self.id

    def json_view(self):
      return {  'id'          : self.id,
                'username'    : self.username,
                'email'       : self.email,
                'experience'  : self.experience,
                # 'created'     : self.created,

                'feedback_0'  : self.feedback_0,
                'feedback_1'  : self.feedback_1,
                'feedback_2'  : self.feedback_2,
                'feedback_3'  : self.feedback_3,

                # Dont need to return feedback as we don't show the defaults anyway
                'first_run'   : self.first_run,
                'api_key'     : self.api_key }
