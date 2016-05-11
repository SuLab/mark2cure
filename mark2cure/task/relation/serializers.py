import json
from pprint import pprint
from collections import Counter

from .models import Relation, Concept, ConceptText, ConceptDocumentRelationship, RelationAnnotation

from rest_framework import serializers


class RelationSerializer(serializers.ModelSerializer):

    def __init__(self, *args, **kwargs):
        # Don't pass the 'fields' arg up to the superclass
        context = kwargs.pop('context', {})
        user = context.get('user', None)

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
        context = kwargs.pop('context', {})
        user = context.get('user', None)

        super(RelationFeedbackSerializer, self).__init__(*args, **kwargs)

    concepts = serializers.SerializerMethodField()

    def get_concepts(self, relation):

        json_data = {
          "c_1_broken": {
            "id": "zl4RlTGwZM9Ud3CCXpU2VZa7eQVnJj0MdbsRBMGy",
            "text": " is not a "
          },

          "c_2_broken": {
            "id": "RdKIrcaEOnM4DRk25g5jAfeNC6HSpsFZaiIPqZer",
            "text": " is not a "
          },

          "c_d": [{
            "id": "8qota4u8hwtcyp65kz9zm0vjyuxwjt12sko084sn",
            "text": "relates to",
            "desc": "The drug has some type of relation to the disease.",
            "example": "Evidence indicates that Amblyomin-X could be a promising candidate for cancer therapy.",
            "children": [{
                "id": "72zuw4bgzz5dniepb3rmls23nsltwocbk274c98m",
                "text": "exacerbates",
                "desc": "The drug worsens or exacerbates the disease in any manner.",
                "example": "St. John's Wort has been known to worsen symptoms of Parkinson's disease."
              }, {
                "id": "jilhvc5p2cy0atls8659a1fggjvvkmahwuspy2kr",
                "text": "treats",
                "desc": "The drug improves or treats the disease in any manner.",
                "example": "Treatment with gingerol may provide a novel therapeutic strategy for prion-mediated neurotoxicity"
              }, {
                "id": "yrrb92b8vtjmcagjj4nx43sbj8wey2moqagk9ea5",
                "text": "increases risk of",
                "desc": "Administration of the drug increases chances of getting the disease. Considered an 'adverse side effect' of the drug.",
                "example": "Parkinson's disease risk factors include ageing, exposure to toxins and genetic defects."
              }, {
                "id": "jyiczzhhupcp7cmebl422ax5dxe1jkwuq647oaw2",
                "text": "may cause",
                "desc": "Administration of the drug may cause or lead to the disease.",
                "example": "43% of mice given the drug NVP-BHG712 developed pulmonary metastases."
              }, {
                "id": "lt18qfxd1ehj7ymxb29wrv6qa41mocwe6eor9dna",
                "text": "prevents",
                "desc": "Administration of the drug prevents the disease",
                "example": "Malarone is taken daily by travelers wanting to prevent malaria."
              }, {
                "id": "nircjx48im90gy5uzqzc79gk3eagg7m36pti8oqi",
                "text": "other relation, or relation unclear"
              }]
            }, {
            "id": "4mzrh5ub3nla6b1ostx7qdparjl3lrd9o567ubif",
            "text": "has no relation to",
            "desc": "There is no relation between the drug and disease.",
            "example": "To study soft tissue sarcoma, an experimental group of mice was given the drug NVP-BHG712."
          }, {
            "id": "ac1m2h0vxye2vuzv6cljr7be12o0d9oclir0kmr8",
            "text": "cannot be determined"
          }],

          "g_d": [{
            "id": "qq84lkjfh46gmx4a9n1jpwxwrmbajsy1qctb9u8j",
            "text": "relates to",
            "desc": "The gene is related to the disease in some manner.",
            "example": "The gene responsible for Triple A syndrome, AAAS, has recently been identified.",
            "children": [{
              "id": "u0q779rcrevnu6aki694dqka4fnfwvwqgpl06ybl",
              "text": "altered expression is associated with",
              "desc": "Gene expression is the process by which information from a gene is used in the synthesis of a functional gene product (mRNA, protein, or microRNAs). When gene expression is altered, gene expression is either increased or decreased.",
              "example": "Several studies revealed significantly higher EPHB4 expression in malignancies such as prostate cancer."
              }, {
              "id": "04110gzdcxz8niuv83ev08ut7lv0xep4iym5sxm5",
              "text": "mutation is associated with",
              "desc": "Gene mutations or aberrations are permanent changes in the gene DNA sequence, differing from the sequence found in most people. Mutations range in size; they can affect anywhere from a single DNA base pair to a large piece of a chromosome that includes multiple genes.",
              "example": "Mutations in the COL5A or COL3A genes are only a few of the genetic causes of Ehlers-Danlos syndrome."
              }, {
              "id": "rhkmksv5jh0vn7p47uk3fwdior6mlgaubwh1l6ow",
              "text": "and post-translational modifications are associated with",
              "desc": "Translation is protein synthesis. When changes are made to a protein during or after translation, it is considered a post-translational modification.",
              "example": "A small interfering RNA causes knockdown of ATP2C1 expression, resulting in defects in both post-translational processing of wild-type thyroglobulin"
            }, {
              "id": "a86ujdunjj2gt0yeyyl027pp4yqbabotvmu0flwt",
              "text": "other relation, or relation unclear"
            }]}, {
            "id": "7aqs9bklotxhbq3r5dcofvsskiefb1yn2nkt1y4a",
            "text": "has no relation to",
            "desc": "The gene does not relate to the disease.",
            "example": "The precise role of Ngly1 in the ERAD process remains unclear in mammals."
          }, {
            "id": "52d80rv4t0h0g14gb83oamjfm8h9rz19zl1ubzku",
            "text": "cannot be determined"
          }],

          "g_c": [{
            "id": "txh8mu2mrnrffik893gr5h0ir7b1y7plgw94n4j7",
            "text": "relates to",

            "children": [{
              "id": "am1wc2yvdcvwcb3yi298xesplbdktzku6wis49iw",
              "text": "affects"
            }, {
              "id": "5ex6vuro19zeneiwlc8yze6dsq1coxvlpojolwgy",
              "text": "is affected by"
            }, {
              "id": "xdnolju6wacvakqnmz6237zwbh0ta3ftw8mdcp50",
              "text": "other relation, or relation unclear"
            }]

            }, {
            "id": "5mgmlk7rnbbd8q6amuapd3elwjpnbho0raegv59c",
            "text": "has no relation to"
            }, {
              "id": "ytbg2u85xp0hayhv1sh6kjlonaag1n4kfgidp1o1",
              "text": "cannot be determined"
            }
          ]
        }
        # relation = 152207 & document = 2621 testing

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

        def get_label_and_group_from_json(json_data):

                answer_tally = Counter(relation_answer_list)
                parent_series_dict_list = []
                children_series_dict_list = []
                for key, value in json_data.items():
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
                            value = [answer_tally[dict_item['id']]]
                        except:
                            value = [0]
                        parent_series_dict_list.append({'label': label, 'values': value, 'group': group})

                        children = check_for_children(dict_item)
                        if children is not None:
                            for dict_in_children in children:
                                label = dict_in_children['text']
                                # group = 'Relates to'
                                group = ''
                                try:
                                    value = [answer_tally[dict_in_children['id']]]
                                except:
                                    value = [0]
                                children_series_dict_list.append({'label': label, 'values': value, 'group': group})

                series_dict_list = parent_series_dict_list + children_series_dict_list
                return series_dict_list

        series_dict_list = get_label_and_group_from_json(json_data)

        return {
            'series': series_dict_list,
            'concept_1_text': cdr1.concept_text.text,
            'concept_2_text': cdr2.concept_text.text
        }

    class Meta:
        model = Relation
