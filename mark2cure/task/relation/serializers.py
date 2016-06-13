from collections import Counter

from .models import Relation, ConceptDocumentRelationship, RelationAnnotation
from . import relation_data

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


class RelationFeedbackSerializer(serializers.ModelSerializer):

    def __init__(self, *args, **kwargs):
        super(RelationFeedbackSerializer, self).__init__(*args, **kwargs)

    concepts = serializers.SerializerMethodField()

    def get_concepts(self, relation):

        def check_for_children(dict_item):
            try:
                return dict_item['children']
            except KeyError:
                return

        def create_list(value):
            if type(value) is list:
                return value
            if type(value) is dict:
                return [value]

        relation_answers = RelationAnnotation.objects.filter(relation=relation)

        relation_answer_list = []
        for i in relation_answers:
            relation_answer_list.append(i.answer)

        cdr_query = ConceptDocumentRelationship.objects.filter(document=relation.document)
        # (TODO) check that this doesn't contain more concept texts

        cdr1 = cdr_query.filter(concept_text__concept_id=relation.concept_1).first()
        cdr2 = cdr_query.filter(concept_text__concept_id=relation.concept_2).first()

        def get_label_and_group_from_json(relation_data):

                answer_tally = Counter(relation_answer_list)
                parent_series_dict_list = []
                children_series_dict_list = []
                for key, value in relation_data.items():
                    dict_list = create_list(value)

                    for dict_item in dict_list:
                        label = dict_item['text']
                        # group = 'Relation'
                        group = ''
                        if key == 'c_1_broken':
                            label = cdr1.concept_text.text + ' marked incorrectly'
                            pass
                        elif key == 'c_2_broken':
                            label = cdr2.concept_text.text + ' marked incorrectly'
                            pass
                        elif key != relation.relation_type:
                            continue

                        try:
                            value = answer_tally[dict_item['id']]
                        except:
                            value = 0
                        parent_series_dict_list.append({'label': label, 'value': value, 'group': group})

                        children = check_for_children(dict_item)
                        if children is not None:
                            for dict_in_children in children:
                                label = dict_in_children['text']
                                # group = 'Relates to'
                                group = ''
                                try:
                                    value = answer_tally[dict_in_children['id']]
                                except:
                                    value = 0
                                children_series_dict_list.append({'label': label, 'value': value, 'group': group})

                series_dict_list = parent_series_dict_list + children_series_dict_list
                return series_dict_list

        series_dict_list = get_label_and_group_from_json(relation_data)

        return {
            'data': series_dict_list,
            'concept_1_text': cdr1.concept_text.text,
            'concept_2_text': cdr2.concept_text.text
        }

    class Meta:
        model = Relation
