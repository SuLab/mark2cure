from django.conf.urls import patterns, include, url
from django.contrib import admin
from django.contrib.flatpages import views

urlpatterns = patterns(
    '',
    url(r'', include('mark2cure.common.urls')),

    url(r'^account/', include('mark2cure.account.urls', namespace='account', app_name='account')),
    url(r'^document/', include('mark2cure.document.urls')),

    url('', include('social.apps.django_app.urls', namespace='social')),
    url('', include('django.contrib.auth.urls', namespace='auth')),

    url(r'^api-auth/',
        include('rest_framework.urls', namespace='rest_framework')),
    url(r'^grappelli/', include('grappelli.urls')),
    url(r'^admin/', include(admin.site.urls)),

    url(r'^(?P<url>.*/)$', views.flatpage),
)
