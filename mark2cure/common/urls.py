from django.conf import settings
from django.conf.urls import patterns, include, url


urlpatterns = patterns(
    'mark2cure.common.views',
    url(r'^$', r'home'),

    # Initial training for fresh signups
    url(r'^training/basics/$', r'introduction'),
    url(r'^training/intro/1/step/(?P<step_num>\w+)/$',
        r'training_one'),
    url(r'^training/intro/2/step/(?P<step_num>\w+)/$',
        r'training_two'),
    url(r'^training/intro/$', r'training_read'),


    url(r'^dashboard/$', r'dashboard'),

    url(r'^message/$', r'message'),
    url(r'^profile-survey/$', r'profile_survey'),
    url(r'^survey/$', r'survey')
)
