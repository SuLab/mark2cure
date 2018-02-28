# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('relation', '0005_auto_20151215_1549'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='relation',
            name='concept_text_1',
        ),
        migrations.RemoveField(
            model_name='relation',
            name='concept_text_2',
        ),
        migrations.AddField(
            model_name='relation',
            name='concept_1',
            field=models.ForeignKey(related_name='concept_1', default=1, to='relation.Concept'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='relation',
            name='concept_2',
            field=models.ForeignKey(related_name='concept_2', default=1, to='relation.Concept'),
            preserve_default=False,
        ),
    ]
