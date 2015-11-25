# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import picklefield.fields


class Migration(migrations.Migration):

    dependencies = [
        ('common', '__first__'),
    ]

    operations = [
        migrations.CreateModel(
            name='Report',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('report_type', models.CharField(max_length=1, choices=[(0, b'pairwise'), (1, b'average')])),
                ('dataframe', picklefield.fields.PickledObjectField(editable=False)),
                ('args', picklefield.fields.PickledObjectField(editable=False)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('group', models.ForeignKey(blank=True, to='common.Group', null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
