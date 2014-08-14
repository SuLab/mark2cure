from django.conf.urls import patterns, include, url


urlpatterns = patterns('',
    (r'^logout/$', 'django.contrib.auth.views.logout', {'next_page': '/'}),
    (r'^login/$', 'django.contrib.auth.views.login', {'template_name': 'account/login.jade'}),
    (r'^reset/$', 'django.contrib.auth.views.password_reset', {'template_name': 'account/reset.jade', 'post_reset_redirect': '/account/reset-thanks/'}),
    (r'^confirm/(?P<uidb36>[0-9A-Za-z]+)-(?P<token>.+)/$', 'django.contrib.auth.views.password_reset_confirm', {'template_name': 'account/confirm.jade', 'post_reset_redirect': '/account/login/'}),
    (r'^change/$', 'django.contrib.auth.views.password_change', {'template_name': 'account/change.jade', 'post_change_redirect': '/account/settings/'}),
)

urlpatterns += patterns('mark2cure.account.views',
    url(r'^reset-thanks/$', r'reset_thanks', name='reset_thanks'),
    url(r'^settings/$', r'settings'),
    url(r'^create/$', r'create'),
    url(r'^$', r'settings'),
    url(r'^(?P<profile_id>\d+)/profile/$', r'update_profile'),
    url(r'^(?P<user_id>\d+)/inactivate/$', r'inactivate'),
    url(r'^(?P<user_id>\d+)/activate/$', r'activate'),
    url(r'^(?P<user_id>\d+)/delete/$', r'delete'),
    url(r'^(?P<user_id>\d+)/staffify/$', r'staffify'),
    url(r'^(?P<user_id>\d+)/destaffify/$', r'destaffify'),
)
