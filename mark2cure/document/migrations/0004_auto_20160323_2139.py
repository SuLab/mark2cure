# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('document', '0003_auto_20151222_1116'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='annotation',
            name='start',
        ),
        migrations.RemoveField(
            model_name='annotation',
            name='text',
        ),
        migrations.RemoveField(
            model_name='annotation',
            name='type',
        ),
        migrations.AlterField(
            model_name='annotation',
            name='kind',
            field=models.CharField(default=b'e', max_length=1, choices=[(b'e', b'Entity Recognition'), (b'r', b'Relation')]),
            preserve_default=True,
        ),
    ]
