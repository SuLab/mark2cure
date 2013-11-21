'''
URLs for static pages
e.g. an about, FAQ, or help page.
'''

from django.conf import settings
from django.conf.urls import patterns, include, url

urlpatterns = patterns('mark2cure.common.views',
    url(r'^$', r'home'),
)
