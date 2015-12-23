# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('relation', '0004_auto_20151215_1228'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='relation',
            name='concept1_id',
        ),
        migrations.RemoveField(
            model_name='relation',
            name='concept2_id',
        ),
        migrations.RemoveField(
            model_name='relation',
            name='relation',
        ),
        migrations.AddField(
            model_name='relation',
            name='concept_text_1',
            field=models.ForeignKey(related_name='concept_text_1', default=1, to='relation.Concept'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='relation',
            name='concept_text_2',
            field=models.ForeignKey(related_name='concept_text_2', default=1, to='relation.Concept'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='relation',
            name='relation_type',
            field=models.CharField(default='', max_length=3),
            preserve_default=False,
        ),
    ]
