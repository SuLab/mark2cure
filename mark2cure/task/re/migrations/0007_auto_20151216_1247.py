# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('relation', '0006_auto_20151215_1821'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='answer',
            name='relation',
        ),
        migrations.DeleteModel(
            name='Answer',
        ),
    ]
