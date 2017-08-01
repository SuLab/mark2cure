from rest_framework import serializers
from . import relation_data_flat


class DocumentRelationSerializer(serializers.Serializer):
    """Organize the flat SQL response into a slightly
        nested response for a Document's available
        relation tasks for a specific user
    """

    id = serializers.IntegerField()
    document_id = serializers.IntegerField()
    relation_type = serializers.CharField()

    community_progress = serializers.DecimalField(max_digits=3, decimal_places=2, coerce_to_string=False)
    community_completed = serializers.BooleanField()
    user_completed = serializers.BooleanField()

    concepts = serializers.SerializerMethodField()

    def get_concepts(self, relation):
        return {
            'c1': {
                'text': relation.get('concept_1_text'),
                'type': relation.get('concept_1_type'),
                'id': relation.get('concept_1_id')
            },
            'c2': {
                'text': relation.get('concept_2_text'),
                'type': relation.get('concept_2_type'),
                'id': relation.get('concept_2_id')
            }
        }


class RelationAnswerSerializer(serializers.Serializer):
    user_id = serializers.IntegerField()
    self = serializers.BooleanField()

    answer = serializers.SerializerMethodField()

    def get_answer(self, obj):
        return next(d for d in relation_data_flat if d['id'] == obj.get('answer_hash'))


class RelationAnalysisSerializer(serializers.Serializer):
    id = serializers.IntegerField()
    document_id = serializers.IntegerField()
    kind = serializers.CharField()

    concept_a = serializers.SerializerMethodField()
    concept_b = serializers.SerializerMethodField()

    answers = RelationAnswerSerializer(many=True)

    def get_concept_a(self, obj):
        return {
            'id': obj.get('concept_1_id'),
            'text': obj.get('concept_1_text')
        }

    def get_concept_b(self, obj):
        return {
            'id': obj.get('concept_2_id'),
            'text': obj.get('concept_2_text')
        }


