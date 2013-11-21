from django.db import models
from django.contrib.auth.models import User

from timezone_field import TimeZoneField
import datetime

class UserProfile(models.Model):
    user                  = models.OneToOneField(User)
    created_by            = models.ForeignKey(User, null=True, blank=True, related_name="children")
    timezone              = TimeZoneField(default='America/Los_Angeles')
    instructions_enabled  = models.BooleanField(default=True, verbose_name="Display Extra Instructions")

    experience  = models.IntegerField(default=0)
    feedback_0  = models.IntegerField(default=0)
    feedback_1  = models.IntegerField(default=0)
    feedback_2  = models.IntegerField(default=0)
    feedback_3  = models.IntegerField(default=0)

    first_run   = models.BooleanField(default = False, blank = True)
    email_bool  = models.BooleanField(default = False, blank = True)

    mturk     = models.BooleanField(default = False, blank = True)

#     views       = db.relationship('View',         backref=db.backref('user',  lazy='select'))
#     ncbos       = db.relationship('Ncbo',         backref=db.backref('user',  lazy='select'))
#     messages    = db.relationship('Message',      backref=db.backref('user',  lazy='select'))
#     annotations = db.relationship('Annotation',   backref=db.backref('user',  lazy='select'))
#     quests      = db.relationship('Quest',        backref=db.backref('user',  lazy='select'))

class Message(models.Model):
    message = models.TextField(blank=True)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    # user_id = db.Column(db.Integer, db.ForeignKey('user.id'))

class Ncbo(models.Model):
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    min_term_size = models.IntegerField()
    score         = models.IntegerField()

    # user_id    = db.Column(db.Integer, db.ForeignKey('user.id'))

    # def __init__(self, min_term_size, score):
      # username = 'annotator_bot_' + str(min_term_size) + '_'+ str(score)

      # user = db.session.query(User).filter_by(username = username).first()
      # if user is None:
      #   user = User(username,
      #               None,
      #               None,
      #               False,
      #               None)
      #   db.session.add(user)
      #   db.session.commit()

      # self.min_term_size = min_term_size
      # self.score = score
      # self.user = user

# class User(db.Model):

#     # Flask-Login integration
#     def is_authenticated(self):
#         return True
#
#     def is_active(self):
#         return True
#
#     def is_anonymous(self):
#         return False
#
#     def get_id(self):
#         return self.id
#
#     def set_email(self, email):
#         self.email = email
#
#     # Required for administrative interface
#     def __unicode__(self):
#         return self.email
#
#     def __repr__(self):
#         return '<User %r>' % self.id
#
#     def __init__(self, username, email=None, experience=0, email_bool=False, mturk=False):
#         self.username   = username
#         self.email      = email
#         self.experience = experience
#         self.created    = datetime.datetime.utcnow()
#
#         self.feedback_0 = -1
#         self.feedback_1 = -1
#         self.feedback_2 = -1
#         self.feedback_3 = -1
#
#         self.first_run  = True
#         self.email_bool = email_bool
#         self.mturk = mturk
#
#     def json_view(self):
#       return {  'id'          : self.id,
#                 'username'    : self.username,
#                 'email'       : self.email,
#                 'experience'  : self.experience,
#                 # 'created'     : self.created,
#
#                 'feedback_0'  : self.feedback_0,
#                 'feedback_1'  : self.feedback_1,
#                 'feedback_2'  : self.feedback_2,
#                 'feedback_3'  : self.feedback_3,
#
#                 # Dont need to return feedback as we don't show the defaults anyway
#                 'mturk'   : self.mturk }
#




User.profile = property(lambda u: UserProfile.objects.get_or_create(user=u)[0])
