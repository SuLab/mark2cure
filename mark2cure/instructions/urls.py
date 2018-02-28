from django.conf.urls import url
from django.views.generic import TemplateView


urlpatterns = [

    # Entity Recognition Task Specific
    url(r'^entity-recognition/disease-marking/$',
        TemplateView.as_view(template_name='instructions/tasks/ner/diseases.html'),
        name='entity-recognition-diseases'),

    url(r'^entity-recognition/gene-marking/$',
        TemplateView.as_view(template_name='instructions/tasks/ner/genes.html'),
        name='entity-recognition-genes'),

    url(r'^entity-recognition/treatment-marking/$',
        TemplateView.as_view(template_name='instructions/tasks/ner/treatments.html'),
        name='entity-recognition-treatments'),

    url(r'^entity-recognition/$',
        TemplateView.as_view(template_name='instructions/tasks/ner/home.html'),
        name='entity-recognition'),

    # Relation Task Specific
    url(r'^relation/definition/disease/$',
        TemplateView.as_view(template_name='instructions/tasks/re/disease-concept.html'),
        name='relation-definition-disease'),

    url(r'^relation/definition/gene/$',
        TemplateView.as_view(template_name='instructions/tasks/re/gene-concept.html'),
        name='relation-definition-gene'),

    url(r'^relation/definition/drug/$',
        TemplateView.as_view(template_name='instructions/tasks/re/drug-concept.html'),
        name='relation-definition-drug'),

    url(r'^relation/gene-disease/$',
        TemplateView.as_view(template_name='instructions/tasks/re/gene-disease-relation.html'),
        name='relation-gene-disease'),

    url(r'^relation/drug-disease/$',
        TemplateView.as_view(template_name='instructions/tasks/re/drug-disease-relation.html'),
        name='relation-drug-disease'),

    url(r'^relation/gene-drug/$',
        TemplateView.as_view(template_name='instructions/tasks/re/gene-drug-relation.html'),
        name='relation-gene-drug'),

    url(r'^relation/$',
        TemplateView.as_view(template_name='instructions/tasks/re/home.html'),
        name='relation'),

    url(r'', TemplateView.as_view(template_name='instructions/home.html'), name='home'),
]
