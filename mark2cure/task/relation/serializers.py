from .models import Relation, ConceptDocumentRelationship

from rest_framework import serializers
from . import relation_data_flat


class RelationCereal(serializers.BaseSerializer):
    def __init__(self, *args, **kwargs):
        context = kwargs.pop('context', {})
        self.sub_dict = context.get('sub_dict', None)
        self.user = context.get('user', None)
        # Instantiate the superclass normally
        super(RelationCereal, self).__init__(*args, **kwargs)

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


class RelationSerializer(serializers.ModelSerializer):

    def __init__(self, *args, **kwargs):
        # Instantiate the superclass normally
        super(RelationSerializer, self).__init__(*args, **kwargs)

    user_completed = serializers.SerializerMethodField('get_user_status')
    concepts = serializers.SerializerMethodField()

    def get_user_status(self, relation):
        return True if relation.user_completed_count > 0 else False

    def get_concepts(self, relation):
        # (TODO) Select the longest text
        cdr_query = ConceptDocumentRelationship.objects.filter(document=relation.document)
        cdr1 = cdr_query.filter(concept_text__concept_id=relation.concept_1).first()
        cdr2 = cdr_query.filter(concept_text__concept_id=relation.concept_2).first()

        return {
            'c1': {
                'text': cdr1.concept_text.text,
                'type': cdr1.stype,
                'id': relation.concept_1.id
            },
            'c2': {
                'text': cdr2.concept_text.text,
                'type': cdr2.stype,
                'id': relation.concept_2.id
            }
        }

    class Meta:
        model = Relation
        fields = ('id', 'document', 'relation_type',
                  'concepts',
                  'user_completed',)

