from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^re/$', views.re_home, name='re'),
    url(r'', views.route, name='route'),
]
