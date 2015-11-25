# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
from django.utils.timezone import utc
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='EmailChangeRequest',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('token', models.CharField(unique=True, max_length=32)),
                ('valid_until', models.DateTimeField(default=datetime.datetime(2015, 12, 2, 21, 6, 32, 793665, tzinfo=utc))),
                ('email', models.EmailField(max_length=75)),
                ('user', models.ForeignKey(related_name='email_change_requests', to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='EmailConfirmationRequest',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('token', models.CharField(unique=True, max_length=32)),
                ('valid_until', models.DateTimeField(default=datetime.datetime(2015, 12, 2, 21, 6, 32, 793665, tzinfo=utc))),
                ('email', models.EmailField(max_length=75)),
            ],
            options={
                'abstract': False,
            },
            bases=(models.Model,),
        ),
    ]
