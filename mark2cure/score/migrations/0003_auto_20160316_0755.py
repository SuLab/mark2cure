# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('score', '0002_auto_20160316_0610'),
    ]

    operations = [
        migrations.AlterField(
            model_name='point',
            name='created',
            field=models.DateTimeField(),
            preserve_default=True,
        ),
    ]
