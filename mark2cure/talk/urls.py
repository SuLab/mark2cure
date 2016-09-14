from django.conf.urls import url
from . import views


urlpatterns = [

    url(r'^(?P<pubmed_id>\d+)/$',
        views.home, name='home'),

    url(r'^annotation/$',
        views.annotation_search, name='annotation-search'),

    url(r'',
        views.recent_discussion, name='recent'),

]
