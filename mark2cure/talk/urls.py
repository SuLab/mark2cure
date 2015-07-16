from django.conf.urls import patterns, include, url
from . import views


urlpatterns = patterns('',


    url(r'^(?P<pubmed_id>\d+)/$',
        views.home, name='home'),
    url(r'^annotation/$',
        views.annotation_search, name='annotation-search'),

    url(r'',
        views.recent_discussion, name='recent'),
)
