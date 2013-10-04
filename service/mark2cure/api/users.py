from flask.ext.restful import reqparse, Resource
from ..core import db
from ..models import User

import os, sys, random, base64, hashlib, string

user_parser = reqparse.RequestParser()
user_parser.add_argument('username',    type=str,   location='json')
user_parser.add_argument('email',       type=str,   location='json')
user_parser.add_argument('experience',  type=int,   location='json')

user_parser.add_argument('feedback_0',  type=int,   location='json')
user_parser.add_argument('feedback_1',  type=int,   location='json')
user_parser.add_argument('feedback_2',  type=int,   location='json')
user_parser.add_argument('feedback_3',  type=int,   location='json')

user_parser.add_argument('first_run',   type=bool,  location='json')
user_parser.add_argument('api_key',     type=str,   location='cookies')

class Users(Resource):
    def get(self, **kwargs):
        args = user_parser.parse_args()
        user = db.session.query(User).filter_by(api_key = args['api_key']).first()
        if not user:
            return {'error' : 'no_user'}, 200
        else:
            # (TODO) return single User object in objects attr?
            return user.json_view(), 200

    def put(self, **kwargs):
        args = user_parser.parse_args()
        user = db.session.query(User).filter_by(api_key = args['api_key']).first()
        if not user:
            return {'error' : 'no_user'}, 200
        else:
            prev_user = db.session.query(User).filter_by(username = args['username']).first()
            if prev_user:
              user = prev_user

            # Set the user keys via the passed Args
            # (TODO) With Args being parsed, can I directly iterate over keys to update?
            user.username   = args['username']
            user.email      = args['email']
            user.experience = args['experience']

            user.feedback_0 = args['feedback_0']
            user.feedback_1 = args['feedback_1']
            user.feedback_2 = args['feedback_2']
            user.feedback_3 = args['feedback_3']

            user.first_run  = args['first_run']

            db.session.commit()
            return user.json_view(), 200

    def post(self, **kwargs):
        # What gets called as soon as we get a new User model
        args = user_parser.parse_args()
        user = User(args['username'],
                    args['email'],
                    args['experience'],
                    base64.b64encode(hashlib.sha256( str(random.getrandbits(256)) ).digest(), random.choice(['rA','aZ','gQ','hH','hG','aR','DD'])).rstrip('=='))
        db.session.add(user)
        db.session.commit()
        return user.json_view(), 201
