from django.db import models
from picklefield.fields import PickledObjectField
from ..common.models import Group


class Report(models.Model):
    # If not group defined, it was ran
    # on the entire annotations as single massive list
    group = models.ForeignKey(Group, null=True, blank=True)

    REPORT_CHOICE = (
        (0, 'pairwise'),
        (1, 'average'),
    )
    report_type = models.CharField(max_length=1, choices=REPORT_CHOICE, blank=False)

    dataframe = PickledObjectField()
    args = PickledObjectField()

    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return u'Report'

