# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='Annotation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('kind', models.CharField(default=b'e', max_length=1, choices=[(b'e', b'Entities'), (b'a', b'Attributes'), (b'r', b'Relations'), (b't', b'Triggers')])),
                ('type', models.CharField(default=b'disease', max_length=40, null=True, blank=True)),
                ('text', models.TextField(null=True, blank=True)),
                ('start', models.IntegerField(null=True, blank=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
            ],
            options={
                'get_latest_by': 'updated',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Document',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('document_id', models.IntegerField(blank=True)),
                ('title', models.TextField()),
                ('authors', models.TextField()),
                ('updated', models.DateTimeField(auto_now=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('source', models.CharField(max_length=200, blank=True)),
            ],
            options={
                'ordering': ('-created',),
                'get_latest_by': 'updated',
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Pubtator',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('kind', models.CharField(max_length=200, blank=True)),
                ('session_id', models.CharField(max_length=200, blank=True)),
                ('content', models.TextField(null=True, blank=True)),
                ('request_count', models.IntegerField(default=0)),
                ('validate_cache', models.BooleanField(default=False)),
                ('updated', models.DateTimeField(auto_now=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('document', models.ForeignKey(to='document.Document')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Section',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('kind', models.CharField(max_length=1, choices=[(b'o', b'Overview'), (b't', b'Title'), (b'a', b'Abstract'), (b'p', b'Paragraph'), (b'f', b'Figure')])),
                ('text', models.TextField(blank=True)),
                ('source', models.ImageField(default=b'images/figure.jpg', upload_to=b'media/images/', blank=True)),
                ('updated', models.DateTimeField(auto_now=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('document', models.ForeignKey(to='document.Document')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='View',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('task_type', models.CharField(default=b'cr', max_length=3, blank=True, choices=[(b'cr', b'Concept Recognition'), (b'cn', b'Concept Normalization'), (b'rv', b'Relationship Verification'), (b'ri', b'Relationship Identification'), (b'rc', b'Relationship Correction')])),
                ('completed', models.BooleanField(default=False)),
                ('opponent', models.ForeignKey(blank=True, to='document.View', null=True)),
                ('section', models.ForeignKey(to='document.Section')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'get_latest_by': 'pk',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='annotation',
            name='view',
            field=models.ForeignKey(to='document.View'),
            preserve_default=True,
        ),
    ]
