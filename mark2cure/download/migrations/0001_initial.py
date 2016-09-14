# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('document', '0004_auto_20160323_2139'),
    ]

    operations = [
        migrations.CreateModel(
            name='Download',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('task_er', models.BooleanField(default=False)),
                ('task_rel', models.BooleanField(default=False)),
                ('file', models.FileField(null=True, upload_to=b'', blank=True)),
                ('create_time', models.DateTimeField(auto_now_add=True)),
                ('download_count', models.IntegerField(default=0)),
                ('documents', models.ManyToManyField(to='document.Document')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
