from django.db import models


class Concept(models.Model):
    document = models.ForeignKey('document.Document')
    uid = models.TextField(blank=False)
    stype = models.TextField(blank=False)
    text = models.TextField(blank=False)

    def __unicode__(self):
        return self.uid


class Relation(models.Model):
    document = models.ForeignKey('document.Document')
    relation = models.TextField(blank=False)
    concept1_id = models.ForeignKey(Concept, related_name='concept1')
    concept2_id = models.ForeignKey(Concept, related_name='concept2')

    def __unicode__(self):
        return self.relation


class Answer(models.Model):
    """
    Class where the user can provide their response to whether or not chemicals
    and diseases are related. TODO: if there is an auto CID, then don't have
    the user provide Answers to the relations!!
    """
    relation = models.ForeignKey(Relation)
    # TODO, this might get updated from Toby's code in the future
    relation_type = models.TextField(blank=True)
    # user confidence order 1 to 4 (where 1 is not confident and 4 is confident)
    # This value is used in scoring later

    username = models.TextField(blank=True)
