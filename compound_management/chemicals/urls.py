
from . import views
from django.urls import path, include

urlpatterns = [
    path('login/', views.login_view, name='login'),
    path('register/', views.register_view, name='register'),
    path('logout/', views.logout_view, name='logout'),
    path('', views.home_view, name='home'),
    path('chemicals/<str:target>/', views.target_view, name='target_view'),
    path('chemicals/<str:target>/new/', views.chemical_new_view, name='chemical_new'),
    path('chemicals/<str:target>/edit/<str:chem_id>/', views.chemical_edit_view, name='chemical_edit'),
    path('chemicals/<str:target>/delete/<str:chem_id>/', views.chemical_delete_view, name='chemical_delete'),
    path('chemicals/upload/<str:target>/', views.upload_chemicals, name='upload_chemicals'),
    path('chemicals/<str:target>/pharmacokinetics/<str:chem_id>/', views.pharmacokinetic_list, name='pharmacokinetic_list'),
    path('chemicals/<str:target>/pharmacokinetics/<str:chem_id>/add/', views.pharmacokinetic_add, name='pharmacokinetic_add'),
    path('chemicals/<str:target>/cytotoxicity/<str:chem_id>/add/', views.cytotoxicity_add, name='cytotoxicity_add'),
    path('chemicals/<str:target>/modeling/<str:chem_id>/', views.schrodinger_model_list, name='schrodinger_model_list'),
    path('chemicals/<str:target>/modeling/<str:chem_id>/add/', views.schrodinger_model_add, name='schrodinger_model_add'),
    path('chemicals/<str:target>/modeling/<str:chem_id>/upload/', views.schrodinger_model_upload, name='schrodinger_model_upload'),
    path('chemicals/<str:target>/liver_stability/<str:chem_id>/add/', views.liver_stability_add, name='liver_stability_add'),
    path('chemicals/<str:target>/cyp_inhibition/<str:chem_id>/add/', views.cyp_inhibition_add, name='cyp_inhibition_add'),
]
