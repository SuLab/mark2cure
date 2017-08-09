from django.conf.urls import url
from . import views


urlpatterns = [

    url(r'^get-started/$', views.get_started, name='get-started'),

    url(r'^dashboard/$',
        views.dashboard, name='dashboard'),

    url(r'^why-i-mark2cure/$',
        views.why_mark2cure, name='why-mark2cure'),

    url(r'^login-with-zooniverse/$',
        views.login_with_zooniverse, name='login-with-zooniverse'),

    url(r'^zooniverse-callback/$',
        views.zooniverse_callback, name='zooniverse-callback'),

    url(r'^group/(?P<group_stub>\w+)/$',
        views.ner_group_home, name='ner-group-home'),

    url(r'^support/$',
        views.support, name='support'),

    url(r'^$',
        views.home, name='home'),

]
