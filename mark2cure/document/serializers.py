from mark2cure.document.models import RelationshipType

from rest_framework import serializers


class RelationshipTypeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = RelationshipType
        fields = ('id', 'full_name', 'type', 'parent')
