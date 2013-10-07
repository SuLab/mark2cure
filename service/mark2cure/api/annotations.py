from flask.ext.restful import reqparse, Resource

class Annotations(Resource):
  def post(self):
    return "post"

  def put(self, ann_id):
    return "put"

  def delete(self, ann_id):
    return "del"
