'''
URLs for static pages
e.g. an about, FAQ, or help page.
'''

from django.conf import settings
from django.conf.urls import patterns, include, url

urlpatterns = patterns('mark2cure.common.views',
    url(r'^$', r'home'),
    url(r'^library/$', r'library'),
    url(r'^library/page/(?P<page_num>\d+)/$', r'library'),
    url(r'^signup/$', r'signup'),
    url(r'^message/$', r'message')
)
