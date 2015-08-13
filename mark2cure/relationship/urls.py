from django.conf.urls import patterns, url

from . import views

# TODO, these are not used currently

urlpatterns = patterns('',
    url(r'^(?P<paper_pk>\d+)/$',
        views.verify_relationship, name='verify_relationship'),

    url(r'^(?P<relation_id>[0-9]+)/relationship_type/$',
        views.relationship_type, name='relationship_type'),

    url(r'',
        views.home, name='home')
)
