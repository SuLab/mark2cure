# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('task', '0003_auto_20160419_2323'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='level',
            options={'ordering': ('-level',), 'get_latest_by': 'created'},
        ),
    ]
