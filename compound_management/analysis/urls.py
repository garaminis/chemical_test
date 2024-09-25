from django.urls import path
from . import views

urlpatterns = [
    path('r/', views.patient_input, name='patient_input' ),
    path('sl/', views.SLselected_gene_input, name='SLselected_gene_input'),
    path('load/', views.r_script_process, name='r_script_process'),
]
