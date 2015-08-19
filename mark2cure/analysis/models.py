from django.db import models
from picklefield.fields import PickledObjectField
from ..common.models import Group


class Report(models.Model):
    # If not group defined, it was ran
    # on the entire annotations as single massive list
    group = models.ForeignKey(Group, null=True, blank=True)

    PAIRWISE = 0
    AVERAGE = 1
    REPORT_CHOICE = (
        (PAIRWISE, 'pairwise'),
        (AVERAGE, 'average'),
    )
    report_type = models.CharField(max_length=1, choices=REPORT_CHOICE, blank=False)

    dataframe = PickledObjectField()
    args = PickledObjectField()

    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return u'Report'

