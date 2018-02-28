# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('relation', '0007_auto_20151216_1247'),
    ]

    operations = [
        migrations.CreateModel(
            name='RelationAnnotation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('answer', models.CharField(max_length=40)),
                ('relation', models.ForeignKey(to='relation.Relation')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
