from django.conf.urls import url
from django.views.generic import TemplateView

from . import views

app_name = 'grna'
urlpatterns = [
    # The initial view for designing new gRNA
    url(r'^$', views.GuideRNAView.as_view()),
    # The results of the query
    url(r'^results/$', views.grna_results, name='grna_results'),
]
