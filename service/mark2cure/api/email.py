from flask.ext.restful import reqparse, Resource
from models import db, User

import os, sys, random, base64, hashlib, string

emailsub_parser = reqparse.RequestParser()
emailsub_parser.add_argument('email',       type=str,   location='form')
emailsub_parser.add_argument('email_bool',  type=str,  location='form')

class Email(Resource):

    def post(self, **kwargs):
        # What gets called as soon as we get a new User model
        args = emailsub_parser.parse_args()

        beta = False
        if args['email_bool'] == "true":
          beta = True


        user = User(args['email'],
                    args['email'],
                    0,
                    base64.b64encode(hashlib.sha256( str(random.getrandbits(256)) ).digest(), random.choice(['rA','aZ','gQ','hH','hG','aR','DD'])).rstrip('=='),
                    beta
                    )
        db.session.add(user)
        db.session.commit()
        return user.json_view(), 201
