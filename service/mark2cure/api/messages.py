from flask.ext.restful import reqparse, Resource
from flask_login import current_user

from ..core import db
from ..models import User, Message

message_parser = reqparse.RequestParser()
message_parser.add_argument('message',  type=str, location='form')

class Messages(Resource):
    def post(self, **kwargs):
        args = message_parser.parse_args()

        message = Message(args['message'], current_user)
        db.session.add(message)
        db.session.commit()
        return {'status' : 'success'}, 200
