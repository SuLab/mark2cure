# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('score', '0003_auto_20160316_0755'),
    ]

    operations = [
        migrations.AlterField(
            model_name='point',
            name='object_id',
            field=models.IntegerField(null=True, blank=True),
            preserve_default=True,
        ),
    ]
