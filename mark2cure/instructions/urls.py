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

    url(r'^relation-identification/$',
        views.relation_identification, name='relation-identification'),
)
