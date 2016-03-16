# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
from django.utils.timezone import utc


class Migration(migrations.Migration):

    dependencies = [
        ('score', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='point',
            name='completed',
        ),
        migrations.AddField(
            model_name='point',
            name='created',
            field=models.DateTimeField(default=datetime.datetime(2016, 3, 16, 13, 10, 45, 926576, tzinfo=utc), auto_now_add=True),
            preserve_default=False,
        ),
    ]
