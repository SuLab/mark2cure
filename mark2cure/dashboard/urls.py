from django.conf.urls import patterns, url
from . import views


urlpatterns = patterns('',
    url(r'^users/training/$',
        views.users_training, name='users_training'),

    url(r'^pubtator/(?P<pk>\d+)/$',
        views.pubtator_actions, name='pubtator'),


    url(r'^doc/(?P<pk>\d+)/reset-pubtator/$',
        views.document_pubtator_actions, name='document_reset_pubtator'),

    url(r'^doc/(?P<pk>\d+)/$',
        views.document_read, name='document'),


    url(r'^group/(?P<pk>\d+)/$',
        views.group_read, name='group'),

    url(r'^group/$',
        views.group_list, name='groups_home')
)

