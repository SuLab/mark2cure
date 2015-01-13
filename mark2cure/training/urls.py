from django.conf.urls import patterns, url
from . import views

urlpatterns = patterns('',
    # Initial training for fresh signups
    url(r'^training/basics/$', views.introduction, name='introduction'),
    url(r'^training/intro/1/step/(?P<step_num>\w+)/$',
        views.training_one, name='training-1'),
    url(r'^training/intro/2/step/(?P<step_num>\w+)/$',
        views.training_two, name='training-2'),
    url(r'^training/intro/3/$',
        views.training_three, name='training-3'),

    url(r'^training/intro/$',
        views.training_read, name='training-index'),
)

