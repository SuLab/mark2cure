from django.db import models


class Concept(models.Model):
    """Unique concept identifiers from Pubtator files which include
    """
    # This can't be primary_key b/c might have duplicates between versions or ontologies
    id = models.CharField(max_length=200, primary_key=True)

    def __unicode__(self):
        return self.id


class ConceptText(models.Model):
    """A Concept can have multiple textual representations
    """
    concept = models.ForeignKey(Concept)
    text = models.CharField(max_length=200)

    def __unicode__(self):
        return self.text


class ConceptDocumentRelationship(models.Model):
    """Track which textual form is used in a specific Document
    """
    concept_text = models.ForeignKey(ConceptText)
    document = models.ForeignKey('document.Document')
    STYPE_CHOICES = (
        ('d', 'Disease'),
        ('g', 'Gene'),
        ('c', 'Chemical')
    )
    stype = models.CharField(max_length=1, choices=STYPE_CHOICES, default='d')

    def __unicode__(self):
        return self.concept_text.text


class Relation(models.Model):
    """Every relationship between all concepts in the Document
    """
    document = models.ForeignKey('document.Document')

    relation_type = models.CharField(max_length=3)

    concept_1 = models.ForeignKey(Concept, related_name='concept_1')
    concept_2 = models.ForeignKey(Concept, related_name='concept_2')

    def __unicode__(self):
        return '{0} >> {1} ({2})'.format(self.concept_1, self.concept_2, self.relation_type)


class RelationAnnotation(models.Model):
    """Additional fields required for Relation task annotations
    """
    # Only access through Document.Annotation.metadata.RelationAnnotation
    relation = models.ForeignKey(Relation, related_name='annotations')
    answer = models.CharField(max_length=40)

    def __unicode__(self):
        return '{0} for {1}'.format(self.answer, self.relation)


class RelationGroup(models.Model):
    """ A group of related documents that will be presented to
        as Relation Tasks

        Notes:
        Strictly from an organizational level. Groups are not intended to
        limit the number of times individual documents have been annotated
    """

    name = models.CharField(max_length=200)
    stub = models.CharField(max_length=200)
    description = models.TextField(blank=True, null=True)
    order = models.DecimalField(default=0, max_digits=3, decimal_places=3)

    enabled = models.BooleanField(default=False)

    documents = models.ManyToManyField('document.Document')

    def __unicode__(self):
        return '{0} ({1} docs)'.format(self.name, self.documents.count())

