# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('document', '0001_initial'),
        ('relation', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='ConceptDocumentRelationship',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('stype', models.CharField(max_length=5, null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='ConceptText',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('text', models.CharField(max_length=200)),
                ('concept', models.ForeignKey(to='relation.Concept')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='conceptdocumentrelationship',
            name='context_text',
            field=models.ForeignKey(to='relation.ConceptText'),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='conceptdocumentrelationship',
            name='section',
            field=models.ForeignKey(to='document.Section'),
            preserve_default=True,
        ),
        migrations.RemoveField(
            model_name='concept',
            name='document',
        ),
        migrations.RemoveField(
            model_name='concept',
            name='stype',
        ),
        migrations.RemoveField(
            model_name='concept',
            name='text',
        ),
        migrations.RemoveField(
            model_name='concept',
            name='uid',
        ),
        migrations.AlterField(
            model_name='concept',
            name='id',
            field=models.CharField(max_length=200, serialize=False, primary_key=True),
            preserve_default=True,
        ),
    ]
