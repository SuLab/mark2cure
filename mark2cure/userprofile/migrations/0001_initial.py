# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import mark2cure.userprofile.models
import django_countries.fields
from django.conf import settings
# import timezone_field.fields


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='Team',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(help_text='You can create a new team.', max_length=255, verbose_name='Team Name', blank=True)),
                ('description', models.TextField(null=True, blank=True)),
                ('created', models.DateTimeField(auto_now_add=True)),
                ('owner', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
        migrations.CreateModel(
            name='UserProfile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('last_seen', models.DateTimeField(null=True, blank=True)),
                # ('timezone', timezone_field.fields.TimeZoneField(default=b'America/Los_Angeles', null=True, blank=True)),
                ('avatar', models.ImageField(default=b'images/default.jpg', upload_to=mark2cure.userprofile.models._content_file_name, blank=True)),
                ('email_notify', models.BooleanField(default=False)),
                ('gender', models.CharField(default=None, max_length=1, null=True, blank=True, choices=[(b'm', b'Male'), (b'f', b'Female')])),
                ('age', models.IntegerField(default=None, null=True, blank=True)),
                ('occupation', models.CharField(max_length=255, blank=True)),
                ('education', models.IntegerField(default=None, null=True, blank=True, choices=[(0, b'Some elementary'), (1, b'Finished elementary'), (2, b'Some high school'), (3, b'Finished high school'), (4, b'Some community college'), (5, b'Finished community college'), (6, b'Some 4-year college'), (7, b'Finished 4-year college'), (8, b'Some masters program'), (9, b'Finished masters program'), (10, b'Some PhD program'), (11, b'Finished PhD program')])),
                ('science_education', models.IntegerField(default=None, null=True, blank=True, choices=[(0, b'Some elementary'), (1, b'Finished elementary'), (2, b'Some high school'), (3, b'Finished high school'), (4, b'Some community college'), (5, b'Finished community college'), (6, b'Some 4-year college'), (7, b'Finished 4-year college'), (8, b'Some masters program'), (9, b'Finished masters program'), (10, b'Some PhD program'), (11, b'Finished PhD program')])),
                ('country', django_countries.fields.CountryField(blank=True, max_length=2)),
                ('referral', models.TextField(verbose_name='I heard about Mark2Cure from', blank=True)),
                ('motivation', models.TextField(verbose_name='I contribute to Mark2Cure because', blank=True)),
                ('quote', models.TextField(verbose_name='Quote / Signature', blank=True)),
                ('rating_votes', models.PositiveIntegerField(default=0, editable=False, blank=True)),
                ('rating_score', models.IntegerField(default=0, editable=False, blank=True)),
                ('team', models.ForeignKey(blank=True, to='userprofile.Team', null=True)),
                ('user', models.OneToOneField(to=settings.AUTH_USER_MODEL)),
            ],
            options={
            },
            bases=(models.Model,),
        ),
    ]
