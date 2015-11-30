from django.conf.urls import patterns, url

from . import views
from django.contrib.auth import views as reset_views

urlpatterns = patterns('',
    url(r'^login/$', 'django.contrib.auth.views.login',
        {'template_name': 'registration/login.jade'},
        name='login'),
    url(r'^logout/$', 'django.contrib.auth.views.logout',
        {'next_page': '/'},
        name='logout'),

    url(r'^create/settings/$', views.user_creation_settings, name='user_creation_settings'),
    url(r'^create/$', views.user_creation, name='user_creation'),

    #url(r'^request_email_confirmation/$', views.request_email_confirmation,
    #    name='request_email_confirmation'),
    #url(r'^confirm_email/(?P<token>\w+)/$', views.confirm_email,
    #    name='confirm_email'),

    #url(r'^request_email_change/$', views.request_email_change,
    #    name='request_email_change'),
    #url(r'^change_email/(?P<token>\w+)/$', views.change_email,
    #    name='change_email'),

    # Reset Email
    url(r'^password_reset/done/$',
        reset_views.password_reset_done,
        {'template_name': 'password-reset/password_reset_done.jade'},
        name='password_reset_done'),
    url(r'^password_reset/$',
        reset_views.password_reset,
        {'template_name': 'password-reset/password_reset.jade',
        'post_reset_redirect': '/registration/password_reset/done/'},
        name='password-reset'),
)
