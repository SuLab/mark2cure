from django.contrib import admin

from mark2cure.account.models import UserProfile, Message, Ncbo

admin.site.register(UserProfile)
admin.site.register(Message)
admin.site.register(Ncbo)
