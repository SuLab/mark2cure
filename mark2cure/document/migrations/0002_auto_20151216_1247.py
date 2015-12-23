# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('relation', '0007_auto_20151216_1247'),
        ('document', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='annotation',
            name='relation',
            field=models.ForeignKey(blank=True, to='relation.Relation', null=True),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='annotation',
            name='view',
            field=models.ForeignKey(blank=True, to='document.View', null=True),
            preserve_default=True,
        ),
    ]
