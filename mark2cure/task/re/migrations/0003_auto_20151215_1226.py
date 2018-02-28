# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('document', '0001_initial'),
        ('relation', '0002_auto_20151215_1105'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='conceptdocumentrelationship',
            name='section',
        ),
        migrations.AddField(
            model_name='conceptdocumentrelationship',
            name='document',
            field=models.ForeignKey(default=1, to='document.Document'),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='conceptdocumentrelationship',
            name='stype',
            field=models.CharField(max_length=1, null=True, choices=[(b'c', b'MESH'), (b'd', b'MEDIC'), (b'g', b'NCBI Gene')]),
            preserve_default=True,
        ),
    ]
