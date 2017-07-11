from django.conf.urls import url
from . import views


urlpatterns = [
    # BioC File for PMID (basic / no annotations)
    url(r'^(?P<document_pk>\d+)/$',
        views.read_document, name='read-document'),

    # BioC File for Pubtator by PMID
    url(r'^pubtator/(?P<pubmed_id>\d+)/',
        views.read_pubtator_bioc, name='read-pubtator-bioc'),

    # BioC File for Specific Pubtator
    url(r'^pubtator/(?P<pk>\d+)/$',
        views.read_pubtator, name='read-pubtator'),
]
