from django.conf.urls import url
from . import views

urlpatterns = [

    url(r'^$',
        views.home, name='home'),

    url(r'^export/$',
        views.start_export, name='export')

]
