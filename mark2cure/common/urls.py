from django.conf import settings
from django.conf.urls import patterns, include, url


urlpatterns = patterns('mark2cure.common.views',
    url(r'^$', r'home'),
    url(r'^library/$', r'library'),

    url(r'^mturk/$', r'mturk'),
    url(r'^router/$', r'router'),

    url(r'^softblock/$', r'softblock'),

    url(r'^library/page/(?P<page_num>\d+)/$', r'library'),
    url(r'^signup/$', r'signup'),
    url(r'^message/$', r'message'),
    url(r'^survey/$', r'survey')
)
