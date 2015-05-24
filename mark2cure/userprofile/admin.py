from django.contrib import admin
from django.db import models
from .models import UserProfile, Team


class TeamAdmin(admin.ModelAdmin):
    search_fields = ['name', 'owner__username']
    list_display = ['name', 'owner', 'created']

    mymodel = models.ForeignKey(Team)


class UserProfileAdmin(admin.ModelAdmin):
    search_fields = ['user__username', 'referral', 'motivation', 'quote']

    list_display = ['user', 'team',
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
