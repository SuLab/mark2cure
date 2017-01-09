from rest_framework import serializers


class AnnotationListSerializer(serializers.ListSerializer):
    pass


class AnnotationSerializer(serializers.Serializer):
    type_id = serializers.IntegerField()
    text = serializers.CharField()
    start = serializers.IntegerField()
    section_pk = serializers.IntegerField()

    def validate_type_id(self, value):
        if value >= 3:
            raise serializers.ValidationError("Annotation type index is out of range")
        return value

    def validate_text(self, value):
        if not value:
            raise serializers.ValidationError("Annotation text can't be empty")
        return value

    class Meta:
        list_serializer_class = AnnotationListSerializer

    # def create(self, validated_data):
    #     return EntityRecognitionAnnotation.objects.create(**validated_data)
    #     data = EntityRecognitionAnnotation.ANNOTATION_TYPE_CHOICE[int(data)]
    # Must trigger this
    # def before_render(self):
    #     """With all of the annotations submitted successfully, mark the view as complete
    #         before returning
    #     """
    #     return HttpResponse(200)
