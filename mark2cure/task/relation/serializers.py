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


class RelationAnalysisSerializer(serializers.BaseSerializer):
    def __init__(self, *args, **kwargs):
        context = kwargs.pop('context', {})
        self.sub_dict = context.get('sub_dict', None)
        self.user = context.get('user', None)
        # Instantiate the superclass normally
        super(RelationAnalysisSerializer, self).__init__(*args, **kwargs)

    def to_representation(self, obj):

        answers = []
        for x in self.sub_dict[obj[0]]:
            answers.append({
                'user_id': x[5],
                'answer': filter(lambda d: d['id'] == x[3], relation_data_flat)[0],
                'self': x[5] == self.user.pk
            })

        return {
            'id': obj[0],
            'document': obj[1],
            'kind': obj[2],
            'concept_a': {
                'id': obj[3],
                'text': obj[5]
            },
            'concept_b': {
                'id': obj[4],
                'text': obj[6]
            },
            'answers': answers
        }



