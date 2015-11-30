# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0001_initial'),
    ]

    operations = [
    ]

    database_operations = [
        migrations.AlterModelTable('Task', 'task_task'),
        migrations.AlterModelTable('UserQuestRelationship', 'task_userquestrelationship'),
        migrations.AlterModelTable('DocumentQuestRelationship', 'task_documentquestrelationship')
    ]

    state_operations = [
        migrations.DeleteModel('Task'),
        migrations.DeleteModel('UserQuestRelationship'),
        migrations.DeleteModel('DocumentQuestRelationship')
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            database_operations=database_operations,
            state_operations=state_operations
        )
    ]
