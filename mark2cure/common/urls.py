from django.conf.urls import patterns, url
from . import views

urlpatterns = patterns('',

    url(r'^$',
        views.home, name='home'),

    url(r'^beta/$',
        views.beta, name='beta'),

    url(r'^why-i-mark2cure/$',
        views.why_mark2cure, name='why-mark2cure'),

    url(r'^dashboard/$',
        views.dashboard, name='dashboard'),

    # Receive the click, look for pattern in the request, the go to views, views tells you to go to the html
    url(r'^relationships/$',
        views.relationships, name='relationships'),

    url(r'^group/(?P<group_stub>\w+)/$',
        views.group_view, name='group'),

    url(r'^support/$',
        views.support, name='support'),

    # Initial training for fresh signups
    url(r'^quest/(?P<quest_pk>\d+)/(?P<doc_idx>\d+)/results/$',
        views.quest_read_doc_results, name='quest-document-results'),

    # BioC File for Results Page
    url(r'^quest/(?P<quest_pk>\d+)/(?P<doc_pk>\d+)/results/(?P<user_pk>\d+).(?P<format_type>\w+)$',
        views.quest_read_doc_results_bioc, name='quest-document-results-bioc'),

    url(r'^quest/(?P<quest_pk>\d+)/(?P<doc_idx>\d+)/$',
        views.quest_read_doc, name='quest-document'),

    url(r'^quest/(?P<quest_pk>\d+)/(?P<document_pk>\d+)/submit/$',
        views.document_quest_submit, name='doc-quest-submit'),

    # Router to:
    # 1) Redirect to specific doc_idx
    # 2) Redirect to quest feedback page
    url(r'^quest/(?P<quest_pk>\d+)/$',
        views.quest_read, name='quest-home'),

    url(r'^quest/(?P<quest_pk>\d+)/feedback/$',
        views.quest_feedback, name='quest-feedback'),

    url(r'^quest/(?P<quest_pk>\d+)/submit/$',
        views.quest_submit, name='quest-submit'),

)
