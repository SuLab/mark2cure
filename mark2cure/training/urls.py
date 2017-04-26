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
    url(r'^relation/1/$', views.relation_training_one, name='relation-training'),
    url(r'^relation/2/$', views.relation_training_two, name='relation-training-two'),
    url(r'^relation/3/$', views.relation_training_three, name='relation-training-three'),

    url(r'', views.route, name='route'),
]
