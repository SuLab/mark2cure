from django.contrib import admin
from .models import UserProfile, Team


admin.site.register(Team)
admin.site.register(UserProfile)
