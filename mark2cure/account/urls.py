from django.conf.urls import patterns, include, url


urlpatterns = patterns(
    '',
    (r'^logout/$',
        'django.contrib.auth.views.logout',
        {'next_page': '/'}),

    (r'^login/$',
        'django.contrib.auth.views.login',
        {'template_name': 'account/login.jade'}),

    url(r'^password-reset/$',
        'django.contrib.auth.views.password_reset',
        kwargs={'template_name': 'account/reset.jade'}, name='password-reset'),

    url(r'^password_reset/done/$',
        'django.contrib.auth.views.password_reset_done',
        kwargs={'template_name': 'account/reset-thanks.jade'}, name='password_reset_done'),

    url(r'^reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
        'django.contrib.auth.views.password_reset_confirm',
        kwargs={'template_name': 'account/reset-confirm.jade'}, name='password_reset_confirm'),

    url(r'^reset/done/$',
        'django.contrib.auth.views.password_reset_complete',
        kwargs={'template_name': 'account/reset-done.jade'}, name='password_reset_complete'),

)

urlpatterns += patterns('mark2cure.account.views',
    url(r'^$', r'settings'),
    url(r'^settings/$', r'settings'),

    # REST Framework
    url(r'^points/$', r'user_points'),

    url(r'^create/$', r'create'),

    url(r'^newsletter/subscribe/$', r'newsletter_subscribe'),

    url(r'^(?P<profile_id>\d+)/profile/$', r'update_profile'),
)
