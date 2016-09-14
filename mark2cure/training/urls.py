from django.conf.urls import url
from . import views

urlpatterns = [
    # Entity Recognition training
    url(r'^entity-recognition/basics/step/(?P<step_num>\w+)/$',
        views.introduction, name='introduction'),
    url(r'^entity-recognition/intro/1/step/(?P<step_num>\w+)/$',
        views.one, name='one'),
    url(r'^entity-recognition/intro/2/step/(?P<step_num>\w+)/$',
        views.two, name='two'),
    url(r'^entity-recognition/intro/3/step/(?P<step_num>\w+)/$',
        views.three, name='three'),
    url(r'^entity-recognition/intro/4/step/(?P<step_num>\w+)/$',
        views.four, name='four'),

    # Relation training
    url(r'^relation/(?P<part_num>\w+)/step/(?P<step_num>\w+)/$', views.relation_training, name='relation-training'),

    url(r'', views.route, name='route'),
]
