# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('relation', '0003_auto_20151215_1226'),
    ]

    operations = [
        migrations.RenameField(
            model_name='conceptdocumentrelationship',
            old_name='context_text',
            new_name='concept_text',
        ),
    ]
