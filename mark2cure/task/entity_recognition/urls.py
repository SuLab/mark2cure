from django.conf.urls import url

from . import views
from .views import NERDocumentSubmissionView

urlpatterns = [

    url(r'^(?P<task_pk>\d+)/(?P<doc_pk>\d+)/results.json',
        views.ner_quest_document_results, name='ner-quest-document-results'),

    url(r'^quest/(?P<quest_pk>\d+)/(?P<document_pk>\d+)/submit/$',
        NERDocumentSubmissionView.as_view(), name='ner-quest-document-submit'),

    url(r'^quest/(?P<quest_pk>\d+)/submit/$',
        views.ner_quest_submit, name='ner-quest-submit'),

    url(r'^quest/(?P<quest_pk>\d+)/$',
        views.ner_quest, name='ner-quest'),
]
