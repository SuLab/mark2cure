from django.conf.urls import patterns, url

from . import views


urlpatterns = patterns('',
    url(r'^$', views.settings, name='settings-root'),
    url(r'^settings/$', views.settings, name='settings'),
    url(r'^points/$', views.user_points, name='points'),
    url(r'^create/settings/$', views.user_creation_settings, name='user_creation_settings'),
    url(r'^create/$', views.user_creation, name='user_creation'),
)
