from django.db import models
from .managers import EntityRecognitionAnnotationManager
from django.forms.models import model_to_dict


class EntityRecognitionAnnotation(models.Model):
    # Only access through Document.Annotation.metadata.RelationAnnotation

    # Disease, Gene, Protein, et cetera...
    ANNOTATION_TYPE_CHOICE = (
        'disease',
        'gene_protein',
        'drug',
    )
    type = models.CharField(max_length=40, blank=True, null=True, default='disease')

    DISEASE = 0
    GENE = 1
    TREATMENT = 2
    TYPE_CHOICES = (
        (DISEASE, 'Disease'),
        (GENE, 'Gene'),
        (TREATMENT, 'Treatment')
    )
    type_idx = models.IntegerField(choices=TYPE_CHOICES, blank=True, null=True)
    '''
    for e in EntityRecognitionAnnotation.objects.all():
        e.type_idx = ['disease', 'gene_protein', 'drug'].index(e.type)
        e.save()
    '''

    text = models.TextField(blank=True, null=True)

    # (WARNING) Different than BioC
    # This is always the start position relative
    # to the section, not the entire document
    start = models.IntegerField(blank=True, null=True)

    objects = EntityRecognitionAnnotationManager()

    def is_exact_match(self, comparing_annotation):
        required_matching_keys = ['start', 'text', 'type']
        self_d = model_to_dict(self)
        compare_d = model_to_dict(comparing_annotation)
        return all([True if self_d[k] == compare_d[k] else False for k in required_matching_keys])
