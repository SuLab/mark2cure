from django.conf.urls import patterns, url

from . import views


urlpatterns = patterns('',
    url(r'^$',
        views.settings, name='settings'),

    url(r'^points/$',
        views.user_points, name='points'),

    url(r'^alert/$',
        views.alert, name='alert'),

    # Public routes
    url(r'^(?P<username>[\w.@+-]+)/edit/',
        views.settings, name='settings'),

    url(r'^(?P<username>[\w.@+-]+)/',
        views.public_profile, name='public-profile'),

)
