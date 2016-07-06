from django.contrib.contenttypes.fields import GenericForeignKey
from django.conf import settings
from django.db import models

AUTH_USER_MODEL = getattr(settings, 'AUTH_USER_MODEL', 'auth.User')


class Point(models.Model):
    user = models.ForeignKey(AUTH_USER_MODEL, related_name='user')
    via = models.ForeignKey(AUTH_USER_MODEL, null=True, blank=True, related_name='via_user')

    from django.contrib.contenttypes.models import ContentType
    content_type = models.ForeignKey(ContentType)
    object_id = models.IntegerField(blank=True, null=True)
    job = GenericForeignKey('content_type', 'object_id')

    amount = models.FloatField()

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ('-updated',)

    def __unicode__(self):
        return '{0}'.format(self.id)
