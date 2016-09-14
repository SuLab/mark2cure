from django.conf.urls import url

from . import views


urlpatterns = [
    url(r'^(?P<teamname>[\w.@+-]+)/',
        views.home, name='home'),
]
