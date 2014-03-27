from django.conf.urls import patterns, include, url

from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    url(r'', include('mark2cure.common.urls')),
    url(r'^account/', include('mark2cure.account.urls')),
    url(r'^document/', include('mark2cure.document.urls')),
    url(r'^analysis/', include('mark2cure.analysis.urls')),

    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
    url(r'^grappelli/', include('grappelli.urls')),

    url(r'^admin/', include(admin.site.urls))
)
