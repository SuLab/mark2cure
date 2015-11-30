# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('document', '0001_initial'),
        ('common', '0002_auto_20151130_0410'),
    ]

    state_operations = [
        migrations.CreateModel(
            name='DocumentQuestRelationship',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('points', models.IntegerField(default=0, max_length=7, blank=True)),
                ('updated', models.DateTimeField(auto_now=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('document', models.ForeignKey(to='document.Document')),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='Task',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=200)),
                ('kind', models.CharField(default=b'q', max_length=1, choices=[(b't', b'Training'), (b'q', b'Quest')])),
                ('completions', models.IntegerField(default=10, null=True, blank=True)),
                ('points', models.IntegerField(default=0, max_length=6, blank=True)),
                ('experiment', models.IntegerField(null=True, blank=True)),
                ('requires_qualification', models.IntegerField(max_length=6, null=True, blank=True)),
                ('provides_qualification', models.IntegerField(max_length=6, null=True, blank=True)),
                ('meta_url', models.CharField(max_length=200, null=True, blank=True)),
                ('updated', models.DateTimeField(auto_now=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('documents', models.ManyToManyField(to='document.Document', through='task.DocumentQuestRelationship', blank=True)),
                ('group', models.ForeignKey(blank=True, to='common.Group', null=True)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='UserQuestRelationship',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('completed', models.BooleanField(default=False)),
                ('score', models.IntegerField(default=5, max_length=7, blank=True)),
                ('updated', models.DateTimeField(auto_now=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('task', models.ForeignKey(to='task.Task')),
                ('user', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
                ('views', models.ManyToManyField(to='document.View')),
            ],
            options={
                'get_latest_by': 'updated',
            },
            bases=(models.Model,),
        ),
        migrations.AddField(
            model_name='task',
            name='users',
            field=models.ManyToManyField(to=settings.AUTH_USER_MODEL, through='task.UserQuestRelationship', blank=True),
            preserve_default=True,
        ),
        migrations.AddField(
            model_name='documentquestrelationship',
            name='task',
            field=models.ForeignKey(to='task.Task'),
            preserve_default=True,
        ),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(state_operations=state_operations)
    ]
