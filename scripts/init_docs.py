#!/usr/bin/env python

from mark2cure.document.models import Document
from mark2cure.document.utils import get_pubmed_documents

print "get ING _pubmed_documents"
print Document.objects.count()

if Document.objects.count() < 10:
  print "get ING _pubmed_documents"
  get_pubmed_documents()

  # doc = Document.objects.create(username='admin')
  # admin.set_password('djJQsRgcrCW9Hh8BmdwgWbGp7P2KxUccqp59D3RX')
  # admin.is_superuser = True
  # admin.is_staff = True
  # doc.save()

