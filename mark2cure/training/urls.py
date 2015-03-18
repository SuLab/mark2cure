from django.conf.urls import patterns, url
from . import views

urlpatterns = patterns('',
    # Initial training for fresh signups
    url(r'^basics/step/(?P<step_num>\w+)/$',
        views.introduction, name='introduction'),

    url(r'^intro/1/step/(?P<step_num>\w+)/$',
        views.one, name='one'),

    url(r'^intro/2/step/(?P<step_num>\w+)/$',
        views.two, name='two'),

    url(r'^intro/3/step/(?P<step_num>\w+)/$',
        views.three, name='three'),

    url(r'^intro/4/step/(?P<step_num>\w+)/$',
        views.four, name='four'),

    url(r'^intro/$',
        views.read, name='index'),
)

