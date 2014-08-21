from django.conf.urls import patterns, include, url
from django.contrib import admin
from django.contrib.flatpages import views

admin.autodiscover()


urlpatterns = patterns('',
    url(r'', include('mark2cure.common.urls')),

    url(r'^account/', include('mark2cure.account.urls')),
    url(r'^document/', include('mark2cure.document.urls')),
    url(r'^analysis/', include('mark2cure.analysis.urls')),

    url('', include('social.apps.django_app.urls', namespace='social')),
    url('', include('django.contrib.auth.urls', namespace='auth')),

    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
    url(r'^grappelli/', include('grappelli.urls')),
    url(r'^admin/', include(admin.site.urls)),
    url(r'^blog/', include('zinnia.urls')),
    url(r'^comments/', include('django.contrib.comments.urls')),

    url(r'^(?P<url>.*/)$', views.flatpage),
)
