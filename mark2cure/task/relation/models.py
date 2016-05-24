from django.db import models


class Concept(models.Model):
    # This can't be primary_key b/c might have duplicates between versions or ontologies
    id = models.CharField(max_length=200, primary_key=True)

    # ETC Meta
    # origin = models.CharField(null=True, blank=True)
    # version = models.DecimalField(null=True, blank=True)
    # ontology_id = models.DecimalField(null=True, blank=True)

    def __unicode__(self):
        return self.id


class ConceptText(models.Model):
    concept = models.ForeignKey(Concept)
    text = models.CharField(max_length=200)

    def __unicode__(self):
        return self.text


class ConceptDocumentRelationship(models.Model):
    concept_text = models.ForeignKey(ConceptText)
    document = models.ForeignKey('document.Document')
    STYPE_CHOICES = (
        ('c', 'Chemical'),
        ('d', 'Disease'),
        ('g', 'Gene')
    )
    stype = models.CharField(max_length=1, choices=STYPE_CHOICES, null=True, blank=False)

    def __unicode__(self):
        return self.concept_text.text


class Relation(models.Model):
    document = models.ForeignKey('document.Document')

    relation_type = models.CharField(max_length=3)

    concept_1 = models.ForeignKey(Concept, related_name='concept_1')
    concept_2 = models.ForeignKey(Concept, related_name='concept_2')

    def __unicode__(self):
        return self.relation_type


class RelationAnnotation(models.Model):
    # Only access through Document.Annotation.metadata.RelationAnnotation
    relation = models.ForeignKey(Relation)
    answer = models.CharField(max_length=40)


class RelationGroup(models.Model):
    '''
        Strictly from an organizational level. Groups are not intended to
        limit the number of times individual documents have been annotated
    '''

    name = models.CharField(max_length=200)
    stub = models.CharField(max_length=200)
    description = models.TextField(blank=True, null=True)
    order = models.DecimalField(default=0, max_digits=3, decimal_places=3)

    enabled = models.BooleanField(default=False)

    documents = models.ManyToManyField('document.Document')

