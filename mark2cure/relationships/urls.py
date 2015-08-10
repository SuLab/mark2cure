from django.conf.urls import patterns, url

from . import views

# TODO, these are not used currently

urlpatterns = patterns('',
    url(r'^(?P<relationship_pk>\d+)/$',
        views.asdfsadfsdf, name='asdfsadfsdf'),

    url(r'^(?P<pk>[0-9]+)/$', views.DetailView.as_view(), name='detail'),
    url(r'^(?P<pk>[0-9]+)/results/$', views.ResultsView.as_view(), name='results'),

    url(r'',
        views.home, name='home')
)
