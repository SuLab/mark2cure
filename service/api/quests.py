from flask import jsonify
from flask.ext.restful import reqparse, Resource

from models import db, Quest

class Quests(Resource):
  def get(self, quest_id):
    quest = db.session.query(Quest).get(quest_id)
    return jsonify(objects=[relation.json_view() for relation in quest.quest_relations])
