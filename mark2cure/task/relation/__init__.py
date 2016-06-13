
# Load the relation data file for quick access
from django.conf import settings
import json


with open(settings.PROJECT_PATH + '/static/js/tasks/relation-data.json') as data_file:
    relation_data = json.load(data_file)


relation_data_flat = []
relation_data_flat.append(relation_data['c_1_broken'])
relation_data_flat.append(relation_data['c_2_broken'])


def looooop(data):
    new_d = data.copy()
    new_d.pop('children', None)
    relation_data_flat.append(new_d)

    for item in data.get('children', []):
        looooop(item)

for key in [u'g_d', u'c_d', u'g_c']:
    comparison_data = relation_data[key]
    for item in comparison_data:
        looooop(item)

