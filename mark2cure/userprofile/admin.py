from django.contrib import admin
from django.db import models
from .models import UserProfile, Team

from django.contrib.auth.admin import UserAdmin
from django.contrib.auth.models import User

UserAdmin.list_display = ('email', 'first_name', 'last_name',
                          'is_active', 'last_login', 'date_joined',
                          'is_staff', 'is_superuser')

admin.site.unregister(User)
admin.site.register(User, UserAdmin)


class TeamAdmin(admin.ModelAdmin):
    search_fields = ['name', 'owner__username']
    list_display = ['name', 'owner', 'created']

    mymodel = models.ForeignKey(Team)


class UserProfileAdmin(admin.ModelAdmin):
    search_fields = ['user__username', 'referral', 'motivation', 'quote']

    list_display = ['user', 'last_seen', 'team',
            'email_notify', 'gender', 'age',
            'occupation', 'education', 'science_education',
            'country', 'referral', 'motivation',
            'quote']

    readonly_fields = ['user', 'team', 'rating',
            'email_notify', 'gender', 'age',
            'occupation', 'education', 'science_education',
            'country', 'referral', 'motivation',
            'quote']

    mymodel = models.ForeignKey(UserProfile)


admin.site.register(Team, TeamAdmin)
admin.site.register(UserProfile, UserProfileAdmin)
