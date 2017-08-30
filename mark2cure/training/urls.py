from django.conf.urls import url
from . import views

urlpatterns = [
    # Relation training
    url(r'^ner/$', views.ner_home, name='ner'),
    url(r'^re/$', views.re_home, name='re'),

    url(r'', views.route, name='route'),
]
