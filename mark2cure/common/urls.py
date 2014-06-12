from django.conf import settings
from django.conf.urls import patterns, include, url


urlpatterns = patterns('mark2cure.common.views',
    url(r'^$', r'home'),
    url(r'^library/$', r'library'),

    url(r'^mturk/$', r'mturk'),

    url(r'^softblock/$', r'softblock'),
    url(r'^banned/$', r'banned'),

    url(r'^library/page/(?P<page_num>\d+)/$', r'library'),
    url(r'^signup/$', r'signup'),
    url(r'^message/$', r'message'),
    url(r'^profile-survey/$', r'profile_survey'),
    url(r'^survey/$', r'survey')
)
