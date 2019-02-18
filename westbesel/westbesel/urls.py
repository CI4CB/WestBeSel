"""
@author: jibril mammeri | heig-vd
"""
from django.conf.urls import include, url
from prediction import views

from django.contrib import admin
admin.autodiscover()

urlpatterns = [
    url(r'^admin/', include(admin.site.urls)),
    url(r'^$', views.IndexView, name='index'),
    url(r'^user/$', views.UserView, name='user'),
    url(r'^help/$', views.HelpView, name='help'),
    url(r'^logout/$', views.logout_view, name='logout'),
    url(r'^results/(?P<proc_id>.+)/$', views.ResultView, name='results'),
    url(r'^download/(?P<proc_id>.+)/$', views.Download, name='download'),
]
