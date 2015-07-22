from django.conf.urls import patterns, url

from . import views


urlpatterns = patterns('',
    url(r'^(?P<teamname>[\w.@+-]+)/',
        views.home, name='home'),
)
