# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('contenttypes', '0001_initial'),
        ('document', '0002_auto_20151216_1247'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='annotation',
            name='relation',
        ),
        migrations.AddField(
            model_name='annotation',
            name='content_type',
            field=models.ForeignKey(blank=True, to='contenttypes.ContentType', null=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='annotation',
            name='object_id',
            field=models.IntegerField(null=True, blank=True),
            preserve_default=True,
        ),
    ]
