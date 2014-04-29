# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'ConceptRelationship'
        db.create_table(u'document_conceptrelationship', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('concept', self.gf('django.db.models.fields.related.ForeignKey')(related_name='actor', to=orm['document.Concept'])),
            ('relationship', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['document.RelationshipType'], null=True, blank=True)),
            ('target', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='target', null=True, to=orm['document.Concept'])),
            ('annotation', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['document.Annotation'], null=True, blank=True)),
            ('validate', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['document.ConceptRelationship'], null=True, blank=True)),
            ('confidence', self.gf('django.db.models.fields.DecimalField')(default=0, max_digits=11, decimal_places=5)),
            ('updated', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
            ('created', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
        ))
        db.send_create_signal(u'document', ['ConceptRelationship'])

        # Adding model 'RelationshipType'
        db.create_table(u'document_relationshiptype', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('full_name', self.gf('django.db.models.fields.CharField')(max_length=80)),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=80)),
            ('definition', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('parent', self.gf('django.db.models.fields.related.ForeignKey')(blank=True, related_name='children', null=True, to=orm['document.RelationshipType'])),
            ('updated', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
            ('created', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
        ))
        db.send_create_signal(u'document', ['RelationshipType'])

        # Adding model 'Concept'
        db.create_table(u'document_concept', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('concept_id', self.gf('django.db.models.fields.TextField')()),
            ('preferred_name', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('definition', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('semantic_type', self.gf('django.db.models.fields.TextField')(blank=True)),
            ('updated', self.gf('django.db.models.fields.DateTimeField')(auto_now=True, blank=True)),
            ('created', self.gf('django.db.models.fields.DateTimeField')(auto_now_add=True, blank=True)),
        ))
        db.send_create_signal(u'document', ['Concept'])

        # Deleting field 'Annotation.concept'
        db.delete_column(u'document_annotation', 'concept_id')


        # Changing field 'Annotation.text'
        db.alter_column(u'document_annotation', 'text', self.gf('django.db.models.fields.TextField')(null=True))

        # Changing field 'Annotation.start'
        db.alter_column(u'document_annotation', 'start', self.gf('django.db.models.fields.IntegerField')(null=True))

        # Changing field 'Annotation.type'
        db.alter_column(u'document_annotation', 'type', self.gf('django.db.models.fields.CharField')(max_length=40, null=True))
        # Adding M2M table for field concepts on 'Section'
        m2m_table_name = db.shorten_name(u'document_section_concepts')
        db.create_table(m2m_table_name, (
            ('id', models.AutoField(verbose_name='ID', primary_key=True, auto_created=True)),
            ('section', models.ForeignKey(orm[u'document.section'], null=False)),
            ('concept', models.ForeignKey(orm[u'document.concept'], null=False))
        ))
        db.create_unique(m2m_table_name, ['section_id', 'concept_id'])


    def backwards(self, orm):
        # Deleting model 'ConceptRelationship'
        db.delete_table(u'document_conceptrelationship')

        # Deleting model 'RelationshipType'
        db.delete_table(u'document_relationshiptype')

        # Deleting model 'Concept'
        db.delete_table(u'document_concept')

        # Adding field 'Annotation.concept'
        db.add_column(u'document_annotation', 'concept',
                      self.gf('django.db.models.fields.related.ForeignKey')(to=orm['common.Concept'], null=True, blank=True),
                      keep_default=False)


        # Changing field 'Annotation.text'
        db.alter_column(u'document_annotation', 'text', self.gf('django.db.models.fields.TextField')(default=False))

        # Changing field 'Annotation.start'
        db.alter_column(u'document_annotation', 'start', self.gf('django.db.models.fields.IntegerField')(default=False))

        # Changing field 'Annotation.type'
        db.alter_column(u'document_annotation', 'type', self.gf('django.db.models.fields.CharField')(default='', max_length=40))
        # Removing M2M table for field concepts on 'Section'
        db.delete_table(db.shorten_name(u'document_section_concepts'))


    models = {
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
        },
        u'document.annotation': {
            'Meta': {'object_name': 'Annotation'},
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'experiment': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'kind': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'player_ip': ('django.db.models.fields.GenericIPAddressField', [], {'max_length': '39', 'null': 'True', 'blank': 'True'}),
            'start': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            'text': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '40', 'null': 'True', 'blank': 'True'}),
            'updated': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'user_agent': ('django.db.models.fields.CharField', [], {'max_length': '150', 'null': 'True', 'blank': 'True'}),
            'view': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['document.View']"})
        },
        u'document.concept': {
            'Meta': {'object_name': 'Concept'},
            'concept_id': ('django.db.models.fields.TextField', [], {}),
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'definition': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'preferred_name': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'semantic_type': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'updated': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'document.conceptrelationship': {
            'Meta': {'object_name': 'ConceptRelationship'},
            'annotation': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['document.Annotation']", 'null': 'True', 'blank': 'True'}),
            'concept': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'actor'", 'to': u"orm['document.Concept']"}),
            'confidence': ('django.db.models.fields.DecimalField', [], {'default': '0', 'max_digits': '11', 'decimal_places': '5'}),
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'relationship': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['document.RelationshipType']", 'null': 'True', 'blank': 'True'}),
            'target': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'target'", 'null': 'True', 'to': u"orm['document.Concept']"}),
            'updated': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'validate': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['document.ConceptRelationship']", 'null': 'True', 'blank': 'True'})
        },
        u'document.document': {
            'Meta': {'object_name': 'Document'},
            'authors': ('django.db.models.fields.TextField', [], {}),
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'document_id': ('django.db.models.fields.IntegerField', [], {'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'source': ('django.db.models.fields.CharField', [], {'max_length': '200', 'blank': 'True'}),
            'title': ('django.db.models.fields.TextField', [], {}),
            'updated': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'document.relationshiptype': {
            'Meta': {'object_name': 'RelationshipType'},
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'definition': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'full_name': ('django.db.models.fields.CharField', [], {'max_length': '80'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'blank': 'True', 'related_name': "'children'", 'null': 'True', 'to': u"orm['document.RelationshipType']"}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '80'}),
            'updated': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'})
        },
        u'document.section': {
            'Meta': {'object_name': 'Section'},
            'cache': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'concepts': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['document.Concept']", 'null': 'True', 'blank': 'True'}),
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            'document': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['document.Document']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'kind': ('django.db.models.fields.CharField', [], {'max_length': '1'}),
            'source': ('django.db.models.fields.files.ImageField', [], {'default': "'images/figure.jpg'", 'max_length': '100', 'blank': 'True'}),
            'text': ('django.db.models.fields.TextField', [], {'blank': 'True'}),
            'updated': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'validate': ('django.db.models.fields.BooleanField', [], {'default': 'False'})
        },
        u'document.view': {
            'Meta': {'object_name': 'View'},
            'created': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'section': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['document.Section']"}),
            'updated': ('django.db.models.fields.DateTimeField', [], {'auto_now': 'True', 'blank': 'True'}),
            'user': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['auth.User']"})
        }
    }

    complete_apps = ['document']