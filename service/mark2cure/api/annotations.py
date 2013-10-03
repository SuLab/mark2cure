from flask.ext.restful import reqparse, Resource

annotation_parser = reqparse.RequestParser()
annotation_parser.add_argument('api_key',       type=str,   location='cookies')

class Annotations(Resource):
  def post(self):
    args = annotation_parser.parse_args()
    return "post"

  def put(self, ann_id):
    args = annotation_parser.parse_args()
    return "put"

  def delete(self, ann_id):
    args = annotation_parser.parse_args()
    return "del"
