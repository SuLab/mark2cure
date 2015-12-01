# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('document', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Answer',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('relation_type', models.TextField(blank=True)),
                ('username', models.TextField(blank=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Concept',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('uid', models.TextField()),
                ('stype', models.TextField()),
                ('text', models.TextField()),
                ('document', models.ForeignKey(to='document.Document')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Relation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('relation', models.TextField()),
                ('concept1_id', models.ForeignKey(related_name='concept1', to='relation.Concept')),
                ('concept2_id', models.ForeignKey(related_name='concept2', to='relation.Concept')),
                ('document', models.ForeignKey(to='document.Document')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='answer',
            name='relation',
            field=models.ForeignKey(to='relation.Relation'),
            preserve_default=True,
        ),
    ]
