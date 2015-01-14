from django.conf.urls import patterns, url
from . import views

urlpatterns = patterns('',
    # Initial training for fresh signups
    url(r'^basics/$',
        views.introduction, name='introduction'),
    url(r'^intro/1/step/(?P<step_num>\w+)/$',
        views.one, name='1'),
    url(r'^intro/2/step/(?P<step_num>\w+)/$',
        views.two, name='2'),
    url(r'^intro/3/$',
        views.three, name='3'),

    url(r'^intro/$',
        views.read, name='index'),
)

