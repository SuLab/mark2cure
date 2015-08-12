from django.conf.urls import patterns, url
from . import views


urlpatterns = patterns('',
    url(r'^user/training/(?P<format_type>\w+)/$',
        views.user_training, name='user_training'),
    url(r'^user/quest-availability/(?P<format_type>\w+)/$',
        views.user_quest_availability, name='user_quest_availability'),

    url(r'^pubtator/(?P<pk>\d+)/$',
        views.pubtator_actions, name='pubtator'),

    url(r'^document/(?P<pk>\d+)/reset-pubtator/$',
        views.document_pubtator_actions, name='document_reset_pubtator'),
    url(r'^document/(?P<pk>\d+)/$',
        views.document_read, name='document'),

    url(r'^group/create/$',
        views.group_create, name='group_create'),
    url(r'^group/(?P<pk>\d+)/analysis/(?P<format_type>\w+)/$',
        views.group_analysis, name='group_analysis'),
    url(r'^group/(?P<pk>\d+)/$',
        views.group_read, name='group'),
    url(r'^group/$',
        views.group_list, name='group_list'),

    url(r'^$',
        views.home, name='home')
)

