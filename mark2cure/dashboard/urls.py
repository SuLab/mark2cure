from django.conf.urls import patterns, url
from . import views


urlpatterns = patterns('',
    url(r'^users/training/$',
        views.users_training, name='users_training'),

    url(r'^group/(?P<pk>\d+)/$',
        views.group_read, name='group'),

    url(r'^group/$',
        views.group_list, name='groups_home')
)

