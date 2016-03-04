from django.conf.urls import url

from . import views

app_name = 'grna'
urlpatterns = [
    # ex: /grna/
    url(r'^$', views.index, name='index'),
    # ex: /grna/5/
    url(r'^(?P<question_id>[0-9]+)/$', views.detail, name='detail'),
    # ex: /grna/5/results/
    url(r'^(?P<question_id>[0-9]+)/results/$', views.results, name='results'),
    # ex: /grna/5/vote/
    url(r'^(?P<question_id>[0-9]+)/vote/$', views.vote, name='vote'),
]
