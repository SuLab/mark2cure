from ..core import db
from .user import User
import datetime

class Ncbo(db.Model):
    id            = db.Column(db.Integer, primary_key=True)
    min_term_size = db.Column(db.Integer)
    score         = db.Column(db.Integer)

    user_id    = db.Column(db.Integer, db.ForeignKey('user.id'))

    def __init__(self, min_term_size, score):
      username = 'annotator_bot_' + str(min_term_size) + '_'+ str(score)

      user = db.session.query(User).filter_by(username = username).first()
      if user is None:
        user = User(username,
                    None,
                    None,
                    False,
                    None)
        db.session.add(user)
        db.session.commit()

      self.min_term_size = min_term_size
      self.score = score
      self.user = user
