from django.conf.urls import patterns, url

from . import views

# TODO, these are not used currently

urlpatterns = patterns('',
    url(r'^(?P<document_pk>\d+)/$',
        views.relation, name='relation'),

    url(r'^(?P<relation_id>[0-9]+)/results/$',
        views.results, name='results'),

    url(r'^test/results/$',
        views.create_post, name='create_post'),

    url(r'',
        views.home, name='home')


)
