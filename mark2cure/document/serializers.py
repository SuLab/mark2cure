from mark2cure.document.models import *
from rest_framework import serializers

class RelationshipTypeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = RelationshipType
        fields = ('full_name', 'type', 'parent')
