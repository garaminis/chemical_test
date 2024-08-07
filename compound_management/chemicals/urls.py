from django.contrib import admin
from . import views
from django.urls import path
from django.conf import settings
from .views import delete_selected_chems
from django.conf.urls.static import static

urlpatterns = [
    path('', views.home_view, name='home'),
    path('login/', views.login_view, name='login'),
    path('register/', views.register_view, name='register'),
    path('logout/', views.logout_view, name='logout'),

    path('chemicals/<str:target>/', views.target_view, name='target_view'),
    path('chemicals/<str:target>/new/', views.chemical_new_view, name='chemical_new'),
    path('chemicals/upload/<str:target>/', views.upload_chemicals, name='upload_chemicals'),
    path('chemicals/<str:target>/edit/<str:chem_id>/', views.chemical_edit_view, name='chemical_edit'),
    path('chemicals/<str:target>/delete/<str:chem_id>/', views.chemical_delete_view, name='chemical_delete'),
    path('chemicals/<str:target>/delete_selected/', delete_selected_chems, name='delete_selected_chems'),

    path('chemicals/<str:target>/result/<str:chem_id>/', views.pharmacokinetic_list, name='pharmacokinetic_list'),
    path('chemicals/<str:target>/result/<str:chem_id>/add/', views.pharmacokinetic_add, name='pharmacokinetic_add'),
    path('chemicals/<str:target>/cytotoxicity/<str:chem_id>/add/', views.cytotoxicity_add, name='cytotoxicity_add'),
    path('chemicals/<str:target>/liver_stability/<str:chem_id>/add/', views.liver_stability_add, name='liver_stability_add'),
    path('chemicals/<str:target>/cyp_inhibition/<str:chem_id>/add/', views.cyp_inhibition_add, name='cyp_inhibition_add'),

    path('chemicals/<str:target>/pharmacokinetic/<str:chem_id>/del/<int:id>', views.pharmacokinetic_delete, name='pharmacokinetic_delete'),
    path('chemicals/<str:target>/cytotoxicity/<str:chem_id>/del/<int:id>', views.cytotoxicity_delete, name='cytotoxicity_delete'),
    path('chemicals/<str:target>/liver_stability/<str:chem_id>/del/<int:id>', views.liver_stability_delete, name='liver_stability_delete'),
path('chemicals/<str:target>/CYPInhibition/<str:chem_id>/del/<int:id>', views.cyp_inhibition_delete, name='cyp_inhibition_delete'),

    path('chemicals/<str:target>/modeling/<str:chem_id>/', views.schrodinger_model_list, name='schrodinger_model_list'),
    path('chemicals/<str:target>/modeling/<str:chem_id>/add/', views.schrodinger_model_add, name='schrodinger_model_add'),
    path('chemicals/<str:target>/modeling/<str:chem_id>/upload/', views.schrodinger_model_upload, name='schrodinger_model_upload'),

    path('r/', views.patient_input, name='patient_input'),
    path('sl/', views.SLselected_gene_input, name='SLselected_gene_input'),

    path('chemicals/<str:target>/cck/<str:chem_id>/add/',views.cck_add, name='cck_add'),
    path('chemicals/<str:target>/wb/<str:chem_id>/add/',views.wb_add, name='wb_add'),
    path('chemicals/<str:target>/in_target/<str:chem_id>/add/',views.in_target_add, name='in_target_add'),
    path('chemicals/<str:target>/other/<str:chem_id>/add/',views.other_add, name='other_add'),

    path('chemicals/<str:target>/cck/<str:chem_id>/del/<int:id>', views.cck_delete, name='cck_delete'),
    path('chemicals/<str:target>/wb/<str:chem_id>/del/<int:id>', views.wb_delete, name='wb_delete'),
    path('chemicals/<str:target>/it/<str:chem_id>/del/<int:id>', views.in_target_delete, name='in_target_delete'),
    path('chemicals/<str:target>/ot/<str:chem_id>/del/<int:id>', views.other_delete, name='other_delete'),
]
if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
