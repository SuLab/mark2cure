from django.conf.urls import patterns, include, url
from django.contrib import admin
from django.contrib.flatpages import views

urlpatterns = patterns('',

    url(r'', include('mark2cure.common.urls',
        namespace='common')),

    url(r'^training/', include('mark2cure.training.urls',
        namespace='training')),
    url(r'^document/', include('mark2cure.document.urls',
        namespace='document')),

    url(r'u/', include('mark2cure.userprofile.urls',
        namespace='profile', app_name='userprofile')),
    url(r'^registration/', include('mark2cure.registration.urls',
        namespace='registration', app_name='registration')),

    url('', include('social.apps.django_app.urls', namespace='social')),
    url('', include('django.contrib.auth.urls', namespace='auth')),

    url(r'^api-auth/',
        include('rest_framework.urls', namespace='rest_framework')),
    url(r'^grappelli/', include('grappelli.urls')),
    url(r'^admin/', include(admin.site.urls)),

    url(r'^(?P<url>.*/)$', views.flatpage),
)
