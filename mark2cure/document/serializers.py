from django.contrib.auth.models import User
from mark2cure.document.models import Annotation

from rest_framework import serializers


class AnnotationSerializer(serializers.ModelSerializer):
    class Meta:
        model = Annotation
        fields = ('text', 'start',)


class TopUserFromViewsSerializer(serializers.ModelSerializer):
    id = serializers.RelatedField(source='user')
    annotations = serializers.SerializerMethodField('get_annotations')

    def get_annotations(self, user):
        # cars_queryset = Car.objects.all().filter(Q(garage=garage) | ...).select_related()
        # serializer = CarSerializer(instance=cars_queryset, many=True, context=self.context)
        user_id = user.get('user')

        queryset = Annotation.objects.filter(view__user__id=user_id).all()[:5]
        serializer = AnnotationSerializer(instance=queryset, many=True, context=self.context)
        return serializer.data

    class Meta:
        model = User
        fields = ('id',)

