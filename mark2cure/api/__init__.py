
# Load the relation data file for quick access
from django.conf import settings
import json

with open(settings.PROJECT_PATH + '/static/js-src/training/training-data.json') as data_file:
    training_data = json.load(data_file)
