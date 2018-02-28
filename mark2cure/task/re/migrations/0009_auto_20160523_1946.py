# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('document', '0004_auto_20160323_2139'),
        ('relation', '0008_relationannotation'),
    ]

    operations = [
        migrations.CreateModel(
            name='RelationGroup',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=200)),
                ('stub', models.CharField(max_length=200)),
                ('description', models.TextField(null=True, blank=True)),
                ('order', models.DecimalField(default=0, max_digits=3, decimal_places=3)),
                ('enabled', models.BooleanField(default=False)),
                ('documents', models.ManyToManyField(to='document.Document')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.AlterField(
            model_name='conceptdocumentrelationship',
            name='stype',
            field=models.CharField(max_length=1, null=True, choices=[(b'c', b'Chemical'), (b'd', b'Disease'), (b'g', b'Gene')]),
            preserve_default=True,
        ),
    ]
