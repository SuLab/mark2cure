# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
import datetime
from django.utils.timezone import utc


class Migration(migrations.Migration):

    dependencies = [
        ('registration', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='emailchangerequest',
            name='valid_until',
            field=models.DateTimeField(default=datetime.datetime(2015, 12, 2, 21, 32, 42, 685479, tzinfo=utc)),
            preserve_default=True,
        ),
        migrations.AlterField(
            model_name='emailconfirmationrequest',
            name='valid_until',
            field=models.DateTimeField(default=datetime.datetime(2015, 12, 2, 21, 32, 42, 685479, tzinfo=utc)),
            preserve_default=True,
        ),
    ]
