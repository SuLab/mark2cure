from django.conf.urls import include, url

from django.contrib import admin
from django.contrib.flatpages import views
from django.contrib.auth import views as reset_views

from django.contrib.sitemaps import Sitemap
from django.contrib.sitemaps.views import sitemap
from django.contrib.flatpages.sitemaps import FlatPageSitemap
from django.urls import reverse

from django.contrib.auth.models import User
from mark2cure.common.models import Group
from mark2cure.userprofile.models import Team


class FlatPages(FlatPageSitemap):
    protocol = 'https'
    priority = 1.0
    changefreq = 'monthly'


class QuestDetails(Sitemap):
    protocol = 'https'
    priority = 1.0
    changefreq = 'weekly'

    def items(self):
        return Group.objects.values_list('stub', flat=True)

    def location(self, item):
        return reverse('common:ner-group-home', kwargs={'group_stub': item})


class UserDetails(Sitemap):
    protocol = 'https'
    priority = .5
    changefreq = 'weekly'

    def items(self):
        return User.objects.values_list('username', flat=True)

    def location(self, item):
        return reverse('profile:public-profile', kwargs={'username': item})


class TeamDetails(Sitemap):
    protocol = 'https'
    priority = .6
    changefreq = 'weekly'

    def items(self):
        return Team.objects.values_list('slug', flat=True)

    def location(self, item):
        return reverse('team:home', kwargs={'slug': item})


sitemaps = {
    'flatpages': FlatPages,
    'quest_detail': QuestDetails,
    'user_detail': UserDetails,
    'team_detail': TeamDetails
}

urlpatterns = [
    url(r'^grappelli/', include('grappelli.urls')),
    url(r'^admin/', include(admin.site.urls)),

    url(r'^sitemap\.xml$', sitemap, {'sitemaps': sitemaps},
        name='django.contrib.sitemaps.views.sitemap'),

    # url(r'^robots\.txt', include('robots.urls')),

    # Password / Account based changes
    url(r'^accounts/', include('mark2cure.userprofile.providers.zooniverse.urls')),
    url(r'^accounts/', include('allauth.urls')),

    # Response / Confirm Changes
    url(r'^reset/done/$',
        reset_views.password_reset_complete,
        {'template_name': 'password-reset/password_reset_complete.html'},
        name='reset'),

    url(r'^reset/confirm/(?P<uidb64>[0-9A-Za-z]+)-(?P<token>.+)/$',
        reset_views.password_reset_confirm,
        {'template_name': 'password-reset/password_reset_confirm.html',
         'post_reset_redirect': '/reset/done/'},
        name='password_reset_confirm'),

    url(r'', include('mark2cure.common.urls',
        namespace='common')),

    url(r'^api/', include('mark2cure.api.urls',
        namespace='api')),

    url(r'^download/', include('mark2cure.download.urls',
        namespace='download')),

    url(r'^control/', include('mark2cure.control.urls',
        namespace='control')),

    url(r'^instructions/', include('mark2cure.instructions.urls',
        namespace='instructions')),
    url(r'^training/', include('mark2cure.training.urls',
        namespace='training')),


    # Task Section
    url(r'^task/ner/', include('mark2cure.task.entity_recognition.urls',
        namespace='task-ner')),

    url(r'^task/re/', include('mark2cure.task.relation.urls',
        namespace='task-re')),

    url(r'^talk/', include('mark2cure.talk.urls',
        namespace='talk')),

    url(r'^team/', include('mark2cure.team.urls',
        namespace='team', app_name='team')),
    url(r'^u/', include('mark2cure.userprofile.urls',
        namespace='profile', app_name='userprofile')),


    url('', include('django.contrib.auth.urls', namespace='auth')),

    url(r'^comments/', include('django_comments.urls')),

    url(r'^api-auth/',
        include('rest_framework.urls', namespace='rest_framework')),

    url(r'^(?P<url>.*/)$', views.flatpage),

]
