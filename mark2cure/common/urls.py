from django.conf.urls import patterns, url
from . import views

urlpatterns = patterns('',

    url(r'^$',
        views.home, name='home'),

    url(r'^beta/$',
        views.beta, name='beta'),

    url(r'^why-i-mark2cure/$',
        views.why_mark2cure, name='why-mark2cure'),

    url(r'^dashboard/$',
        views.dashboard, name='dashboard'),

    url(r'^group/(?P<group_stub>\w+)/network/$',
        views.group_network, name='group-network'),

    url(r'^group/(?P<group_stub>\w+)/$',
        views.group_view, name='group'),

    url(r'^support/$',
        views.support, name='support'),

)
