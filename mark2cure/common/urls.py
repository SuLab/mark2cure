from django.conf import settings
from django.conf.urls import patterns, include, url


urlpatterns = patterns('mark2cure.common.views',
    url(r'^$', r'home'),

    url(r'^training/basics/$', r'introduction'),
    url(r'^training/intro/(?P<quest_num>\d+)/step/(?P<step_num>\w+)/$', r'quest'),
    url(r'^training/intro/$', r'quest_read'),

    url(r'^dashboard/$', r'dashboard'),

    url(r'^signup/$', r'signup'),
    url(r'^message/$', r'message'),
    url(r'^profile-survey/$', r'profile_survey'),
    url(r'^survey/$', r'survey')
)
