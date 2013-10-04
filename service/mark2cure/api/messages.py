from flask.ext.restful import reqparse, Resource
from ..core import db
from ..models import User, Message

message_parser = reqparse.RequestParser()
message_parser.add_argument('message',  type=str, location='form')
message_parser.add_argument('api_key',  type=str, location='cookies')

class Messages(Resource):
    def post(self, **kwargs):
        args = message_parser.parse_args()
        user = db.session.query(User).filter_by(api_key = args['api_key']).first()

        message = Message(args['message'], user)
        db.session.add(message)
        db.session.commit()
        return {'status' : 'success'}, 200
