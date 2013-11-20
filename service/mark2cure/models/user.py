from ..core import db
import datetime

class User(db.Model):
    id          = db.Column(db.Integer, primary_key=True)

    username    = db.Column(db.String(120), unique=True)
    email       = db.Column(db.String(120))
    experience  = db.Column(db.Integer())

    created     = db.Column(db.DateTime)

    feedback_0  = db.Column(db.Integer())
    feedback_1  = db.Column(db.Integer())
    feedback_2  = db.Column(db.Integer())
    feedback_3  = db.Column(db.Integer())

    first_run   = db.Column(db.Boolean())
    email_bool  = db.Column(db.Boolean())

    mturk     = db.Column(db.Boolean())
    admin     = db.Column(db.Boolean())
    active    = db.Column(db.Boolean())

    views       = db.relationship('View',         backref=db.backref('user',  lazy='select'))
    ncbos       = db.relationship('Ncbo',         backref=db.backref('user',  lazy='select'))
    messages    = db.relationship('Message',      backref=db.backref('user',  lazy='select'))
    annotations = db.relationship('Annotation',   backref=db.backref('user',  lazy='select'))
    quests      = db.relationship('Quest',        backref=db.backref('user',  lazy='select'))

    # Flask-Login integration
    def is_authenticated(self):
        return True

    def is_active(self):
        return True

    def is_anonymous(self):
        return False

    def get_id(self):
        return self.id

    def set_email(self, email):
        self.email = email

    # Required for administrative interface
    def __unicode__(self):
        return self.email

    def __repr__(self):
        return '<User %r>' % self.id

    def __init__(self, username, email=None, experience=0, email_bool=False, mturk=False):
        self.username   = username
        self.email      = email
        self.experience = experience
        self.created    = datetime.datetime.utcnow()

        self.feedback_0 = -1
        self.feedback_1 = -1
        self.feedback_2 = -1
        self.feedback_3 = -1

        self.first_run  = True
        self.email_bool = email_bool
        self.mturk = mturk

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
                'mturk'   : self.mturk }
