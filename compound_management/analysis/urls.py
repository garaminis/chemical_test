from django.urls import path
from . import views

urlpatterns = [
    path('analysis/r', views.patient_input, name='patient_input' ),
    path('analysis/sl', views.SLselected_gene_input, name='SLselected_gene_input'),
]
