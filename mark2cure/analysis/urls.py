from django.conf.urls import patterns, url
from . import views


urlpatterns = patterns('',
    url(r'^users/training/$', views.users_training, name='users_training'),
)

