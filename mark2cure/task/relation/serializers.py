from .models import Relation, Concept, ConceptDocumentRelationship

from rest_framework import serializers


class RelationSerializer(serializers.ModelSerializer):

    def __init__(self, *args, **kwargs):
        # Don't pass the 'fields' arg up to the superclass
        context = kwargs.pop('context', {})
        user = context.get('user', None)
        #self.user_highest_level = user.profile.highest_level('skill').level

        # Instantiate the superclass normally
        super(RelationSerializer, self).__init__(*args, **kwargs)

    user_completed = serializers.SerializerMethodField('get_user_status')
    concepts = serializers.SerializerMethodField()

    def get_user_status(self, relation):
        return False

    def get_concepts(self, relation):
        # (TODO) Select the longest text
        cdr_query = ConceptDocumentRelationship.objects.filter(document=relation.document)
        cdr1 = cdr_query.filter(concept_text__concept_id=relation.concept_1).first()
        cdr2 = cdr_query.filter(concept_text__concept_id=relation.concept_2).first()

        return {
            'c1': {
                'text': cdr1.concept_text.text,
                'type': cdr1.stype
            },
            'c2': {
                'text': cdr2.concept_text.text,
                'type': cdr2.stype
            }
        }

    class Meta:
        model = Relation
        fields = ('id', 'document', 'relation_type',
                  'concepts',
                  'user_completed')

