from django.conf.urls import patterns, url
from . import views

urlpatterns = patterns('',
    # Initial training for fresh signups
    url(r'^disease-marking/$',
        views.disease_marking, name='disease-marking'),

    url(r'^gene-marking/$',
        views.gene_marking, name='gene-marking'),

    url(r'^treatment-marking/$',
        views.treatment_marking, name='treatment-marking'),

    url(r'^gene-disease-relation/$',
        views.gene_disease_relation, name='gene-disease-relation'),

    url(r'^chemical-disease-relation/$',
        views.chemical_disease_relation, name='chemical-disease-relation'),

    # TODO: add this after ontology is finalized.
    url(r'^gene-chemical-relation/$',
        views.gene_chemical_relation, name='gene-chemical-relation'),




)
