
# Load the relation data file for quick access
from django.conf import settings
import json


with open(settings.PROJECT_PATH+'/static/js/tasks/relation-data.json') as data_file:
    relation_data = json.load(data_file)
