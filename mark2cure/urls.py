from django.conf.urls import include, url
from django.contrib import admin

from django.contrib.flatpages import views
from django.contrib.auth import views as reset_views

from django.contrib.sitemaps.views import sitemap
from django.contrib.flatpages.sitemaps import FlatPageSitemap


sitemaps = {
    'flatpages': FlatPageSitemap
}

urlpatterns = [
    url(r'^grappelli/', include('grappelli.urls')),
    url(r'^admin/', include(admin.site.urls)),

    url(r'^sitemap\.xml$', sitemap, {'sitemaps': sitemaps},
        name='django.contrib.sitemaps.views.sitemap'),

    # url(r'^robots\.txt', include('robots.urls')),

    # Password / Account based changes
    # url(r'^accounts/signup/$', signup, name='account_signup'),
    url(r'^accounts/', include('mark2cure.userprofile.providers.zooniverse.urls')),
    url(r'^accounts/', include('allauth.urls')),

    # Response / Confirm Changes
    url(r'^reset/done/$',
        reset_views.password_reset_complete,
        {'template_name': 'password-reset/password_reset_complete.jade'},
        name='reset'),

    url(r'^reset/confirm/(?P<uidb64>[0-9A-Za-z]+)-(?P<token>.+)/$',
        reset_views.password_reset_confirm,
        {'template_name': 'password-reset/password_reset_confirm.jade',
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
