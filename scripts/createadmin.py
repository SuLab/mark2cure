#!/usr/bin/env python

from django.contrib.auth.models import User
from mark2cure.settings import ADMIN_PASSWORD

if User.objects.count() == 0:
  admin = User.objects.create(username='admin')
  admin.set_password(ADMIN_PASSWORD)
  admin.is_superuser = True
  admin.is_staff = True
  admin.save()

