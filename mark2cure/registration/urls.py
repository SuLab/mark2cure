from django.conf.urls import patterns, url

from . import views


urlpatterns = patterns('',
    url(r'^login/$', 'django.contrib.auth.views.login',
        {'template_name': 'registration/login.jade',},
        name='login'),
    url(r'^logout/$', 'django.contrib.auth.views.logout',
        {'next_page': '/'},
        name='logout'),

    url(r'^create/settings/$', views.user_creation_settings, name='user_creation_settings'),
    url(r'^create/$', views.user_creation, name='user_creation'),

    url(r'^change_password/$', views.change_password,
        name='change_password'),

    url(r'^request_email_confirmation/$', views.request_email_confirmation,
        name='request_email_confirmation'),
    url(r'^confirm_email/(?P<token>\w+)/$', views.confirm_email,
        name='confirm_email'),

    url(r'^request_email_change/$', views.request_email_change,
        name='request_email_change'),
    url(r'^change_email/(?P<token>\w+)/$', views.change_email,
        name='change_email'),

    # Password reset process
    url(r'^password-reset/$',
        'django.contrib.auth.views.password_reset',
        kwargs={'template_name': 'password-reset/home.jade'},
        name='password-reset'),
    url(r'^password_reset/done/$',
        'django.contrib.auth.views.password_reset_done',
        kwargs={'template_name': 'password-reset/thanks.jade'},
        name='password_reset_done'),
    url(r'^reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$',
        'django.contrib.auth.views.password_reset_confirm',
        kwargs={'template_name': 'password-reset/confirm.jade'},
        name='password_reset_confirm'),
    url(r'^reset/done/$',
        'django.contrib.auth.views.password_reset_complete',
        kwargs={'template_name': 'password-reset/success.jade'},
        name='password_reset_complete'),

    )
