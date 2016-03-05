from django.conf.urls import url

from . import views

app_name = 'grna'
urlpatterns = [
    # The initial view for designing new gRNA
    url(r'^$', views.grna_new, name='grna_new'),
    # The results of the query
    url(r'^results/$', views.grna_results, name='grna_results'),
]
