from django.conf.urls import url
from . import views


urlpatterns = [

    url(r'^$',
        views.home, name='home'),

    url(r'^get-started/$', views.get_started, name='get-started'),

    url(r'^dashboard/$',
        views.dashboard, name='dashboard'),

    url(r'^why-i-mark2cure/$',
        views.why_mark2cure, name='why-mark2cure'),

    # (TODO) this link is broken. Can url routing be removed?
    url(r'^group/(?P<group_stub>\w+)/network/$',
        views.group_network, name='group-network'),

    url(r'^group/(?P<group_stub>\w+)/$',
        views.group_view, name='group'),

    url(r'^support/$',
        views.support, name='support'),

]
