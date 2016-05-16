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

    url(r'^drug-disease-relation/$',
        views.drug_disease_relation, name='drug-disease-relation'),

    url(r'^gene-drug-relation/$',
        views.gene_drug_relation, name='gene-drug-relation'),

)
