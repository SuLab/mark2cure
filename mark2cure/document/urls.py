from django.conf.urls import patterns, url

from . import views


urlpatterns = patterns('',

    # BioC File for PMID (basic / no annotations)
    url(r'^(?P<pubmed_id>\d+).(?P<format_type>\w+)$',
        views.read_bioc, name='read-bioc'),

    # BioC File for Specific Pubtator Response
    url(r'^pubtator/specific/(?P<pub_pk>\d+).(?P<format_type>\w+)$',
        views.read_specific_pubtator_bioc, name='read-specific-pubtator-bioc'),

    # BioC File for Pubtator by PMID
    url(r'^pubtator/(?P<pubmed_id>\d+).(?P<format_type>\w+)$',
        views.read_pubtator_bioc, name='read-pubtator-bioc'),

    url(r'^pubtator-new/(?P<pubmed_id>\d+).json',
        views.read_pubtator_new_bioc, name='read-pubtator-new-bioc'),

    # BioC File for Specific Pubtator
    url(r'^pubtator/(?P<pk>\d+)/$',
        views.read_pubtator, name='read-pubtator'),

)
