# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'UserProfile'
        db.create_table(u'account_userprofile', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['auth.User'], unique=True)),
            ('created_by', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='children', null=True, to=orm['auth.User'])),
            ('timezone', self.gf('timezone_field.fields.TimeZoneField')(default='America/Los_Angeles')),
            ('instructions_enabled', self.gf('django.db.models.fields.BooleanField')(default=True)),
            ('experience', self.gf('django.db.models.fields.IntegerField')(default=0)),
            ('feedback_0', self.gf('django.db.models.fields.IntegerField')(default=0)),
            ('feedback_1', self.gf('django.db.models.fields.IntegerField')(default=0)),
            ('feedback_2', self.gf('django.db.models.fields.IntegerField')(default=0)),
            ('feedback_3', self.gf('django.db.models.fields.IntegerField')(default=0)),
            ('first_run', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('email_notify', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('mturk', self.gf('django.db.models.fields.BooleanField')(default=False)),
            ('turk_last_assignment_id', self.gf('django.db.models.fields.CharField')(max_length=200, blank=True)),
            ('turk_submit_to', self.gf('django.db.models.fields.CharField')(default='http://example.com', max_length=200, blank=True)),
            ('ncbo', self.gf('django.db.models.fields.BooleanField')(default=False)),
        ))
        db.send_create_signal(u'account', ['UserProfile'])

        # Adding model 'Ncbo'
        db.create_table(u'account_ncbo', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('user', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['auth.User'], unique=True, null=True)),
            ('updated', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
            ('created', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
            ('min_term_size', self.gf('django.db.models.fields.IntegerField')()),
            ('score', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'account', ['Ncbo'])


    def backwards(self, orm):
        # Deleting model 'UserProfile'
        db.delete_table(u'account_userprofile')

        # Deleting model 'Ncbo'
        db.delete_table(u'account_ncbo')


    models = {
        u'account.ncbo': {
            'Meta': {'object_name': 'Ncbo'},
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'min_term_size': ('django.db.models.fields.IntegerField', [], {}),
            'score': ('django.db.models.fields.IntegerField', [], {}),
            'updated': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'user': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['auth.User']", 'unique': 'True', 'null': 'True'})
        },
        u'account.userprofile': {
            'Meta': {'object_name': 'UserProfile'},
            'created_by': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'children'", 'null': 'True', 'to': u"orm['auth.User']"}),
            'email_notify': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'experience': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'feedback_0': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'feedback_1': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'feedback_2': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'feedback_3': ('django.db.models.fields.IntegerField', [], {'default': '0'}),
            'first_run': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'instructions_enabled': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'mturk': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'ncbo': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'timezone': ('timezone_field.fields.TimeZoneField', [], {'default': "'America/Los_Angeles'"}),
            'turk_last_assignment_id': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'turk_submit_to': ('django.db.models.fields.CharField', [], {'default': "'http://example.com'", 'max_length': '200', 'blank': 'True'}),
            'user': ('django.db.models.fields.related.OneToOneField', [], {'to': u"orm['auth.User']", 'unique': 'True'})
        },
        u'auth.group': {
            'Meta': {'object_name': 'Group'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '80'}),
            'permissions': ('django.db.models.fields.related.ManyToManyField', [], {'to': u"orm['auth.Permission']", 'symmetrical': 'False', 'blank': 'True'})
        },
        u'auth.permission': {
            'Meta': {'ordering': "(u'content_type__app_label', u'content_type__model', u'codename')", 'unique_together': "((u'content_type', u'codename'),)", 'object_name': 'Permission'},
            'codename': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'content_type': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['contenttypes.ContentType']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'auth.user': {
            'Meta': {'object_name': 'User'},
            'date_joined': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'email': ('django.db.models.fields.EmailField', [], {'max_length': '75', 'blank': 'True'}),
            'first_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'groups': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Group']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_active': ('django.db.models.fields.BooleanField', [], {'default': 'True'}),
            'is_staff': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'is_superuser': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'last_login': ('django.db.models.fields.DateTimeField', [], {'default': 'datetime.datetime.now'}),
            'last_name': ('django.db.models.fields.CharField', [], {'max_length': '30', 'blank': 'True'}),
            'password': ('django.db.models.fields.CharField', [], {'max_length': '128'}),
            'user_permissions': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'related_name': "u'user_set'", 'blank': 'True', 'to': u"orm['auth.Permission']"}),
            'username': ('django.db.models.fields.CharField', [], {'unique': 'True', 'max_length': '30'})
        },
        u'contenttypes.contenttype': {
            'Meta': {'ordering': "('name',)", 'unique_together': "(('app_label', 'model'),)", 'object_name': 'ContentType', 'db_table': "'django_content_type'"},
            'app_label': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'model': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        }
    }

    complete_apps = ['account']