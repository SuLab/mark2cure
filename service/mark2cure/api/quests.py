from flask import jsonify
from flask.ext.restful import reqparse, Resource

from ..core import db
from ..models import Quest

class Quests(Resource):
  def get(self, quest_id=None, quest_name=None):
    if quest_id:
      quest = db.session.query(Quest).get(quest_id)
      return jsonify(objects=[relation.json_view() for relation in quest.quest_relations])
    elif quest_name:
      quest = db.session.query(Quest).filter_by(name = quest_name).first()
      return jsonify(objects=[relation.json_view() for relation in quest.quest_relations])
    else:
      quests = Quest.query.all()
      return jsonify(objects=[quest.json_view() for quest in quests])
