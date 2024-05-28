
from . import views
from django.urls import path, include

urlpatterns = [
    path('login/', views.login_view, name='login'),
    path('register/', views.register_view, name='register'),
    path('logout/', views.logout_view, name='logout'),
    path('', views.home_view, name='home'),
    # path('chemicals/', views.chemical_list_view, name='chemical_list'),
    path('chemicals/<str:target>/', views.target_view, name='target_view'),
    path('chemicals/<str:target>/new/', views.chemical_new_view, name='chemical_new'),
    path('chemicals/<str:target>/edit/<str:chem_id>/', views.chemical_edit_view, name='chemical_edit'),
    path('chemicals/<str:target>/delete/<str:chem_id>/', views.chemical_delete_view, name='chemical_delete'),
    path('chemicals/upload/<str:target>/', views.upload_chemicals, name='upload_chemicals'),
    path('chemicals/<str:target>/pharmacokinetics/<str:chem_id>/', views.pharmacokinetic_list, name='pharmacokinetic_list'),
    path('chemicals/<str:target>/pharmacokinetics/<str:chem_id>/add/', views.pharmacokinetic_add, name='pharmacokinetic_add'),
    path('chemicals/<str:target>/cytotoxicity/<str:chem_id>/add/', views.cytotoxicity_add, name='cytotoxicity_add'),
]
