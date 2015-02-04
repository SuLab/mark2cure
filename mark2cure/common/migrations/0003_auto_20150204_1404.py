# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0002_supportmessage_referral'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='userquestrelationship',
            options={'get_latest_by': 'updated'},
        ),
        migrations.AlterField(
            model_name='supportmessage',
            name='user',
            field=models.ForeignKey(blank=True, to=settings.AUTH_USER_MODEL, null=True),
            preserve_default=True,
        ),
    ]
