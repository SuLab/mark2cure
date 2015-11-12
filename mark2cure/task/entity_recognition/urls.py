from django.conf.urls import patterns, url

from . import views

urlpatterns = patterns('',

    # BioC File for M2C annotations
    url(r'^m2c/(?P<pubmed_id>\d+).(?P<format_type>\w+)$',
        views.read_users_bioc, name='read-users-bioc'),

    # BioC File for Showing User Comparisons accross PMID
    url(r'^(?P<doc_pk>\d+)/user/(?P<user_pk>\d+)/results.(?P<format_type>\w+)$',
        views.user_pmid_results_bioc, name='user-pmid-results-bioc'),

    # BioC File for Partner for Task & PMID Paring
    url(r'^(?P<task_pk>\d+)/(?P<doc_pk>\d+)/results.(?P<format_type>\w+)$',
        views.identify_annotations_results_bioc, name='results-bioc'),

    # "API" Endpoint to submit Document Annotations
    url(r'^(?P<task_pk>\d+)/(?P<section_pk>\d+)/annotation/$',
        views.identify_annotations_submit, name='create'),

)
