from django.contrib import admin
from . import views
from django.urls import path, include
from django.conf import settings
from .views import delete_selected_chems, delete_file
from django.conf.urls.static import static

urlpatterns = [
    path('', views.home_view, name='home'),
# target CRUD
    path('item/<str:chem_id>/toggle_favorite/', views.toggle_favorite, name='toggle_favorite'), # 찜
    path('chemicals/<str:target>/', views.target_view, name='target_view'),
    path('chemicals/<str:target>/new/', views.chemical_new_view, name='chemical_new'),
    path('chemicals/upload/<str:target>/', views.upload_chemicals, name='upload_chemicals'),
    path('chemicals/<str:target>/edit/<str:chem_id>/', views.chemical_edit_view, name='chemical_edit'),
    path('chemicals/<str:target>/delete/<str:chem_id>/', views.chemical_delete_view, name='chemical_delete'),
    path('chemicals/<str:target>/delete_selected/', delete_selected_chems, name='delete_selected_chems'),
# result add,list
    path('chemicals/<str:target>/result/<str:chem_id>/', views.pharmacokinetic_list, name='pharmacokinetic_list'),
    path('chemicals/<str:target>/result/', views.result_add, name='result_add'),
    path('chemicals/<str:target>/result/<str:chem_id>/add/', views.pharmacokinetic_add, name='pharmacokinetic_add'),
    path('chemicals/<str:target>/cytotoxicity/<str:chem_id>/add/', views.cytotoxicity_add, name='cytotoxicity_add'),
    path('chemicals/<str:target>/liver_stability/<str:chem_id>/add/', views.liver_stability_add, name='liver_stability_add'),
    path('chemicals/<str:target>/cyp_inhibition/<str:chem_id>/add/', views.cyp_inhibition_add, name='cyp_inhibition_add'),
    path('chemicals/<str:target>/invivo/<str:chem_id>/add/<int:category>', views.invivo_add, name='invivo_add'),
    path('result/<str:db>/<int:id>/', views.result_file, name='result_file'),
    path('result_img/<str:target>/<str:chem_id>/<int:id>/', views.result_img, name='result_img'),
   # path('chemicals/<str:target>/result2/<str:chem_id>/ajax-pagination/', views.ajax_pagination, name='ajax_pagination'),

    # result delete
    path('chemicals/<str:target>/pharmacokinetic/<str:chem_id>/del/<int:id>', views.pharmacokinetic_delete, name='pharmacokinetic_delete'),
    path('chemicals/<str:target>/cytotoxicity/<str:chem_id>/del/<int:id>', views.cytotoxicity_delete, name='cytotoxicity_delete'),
    path('chemicals/<str:target>/liver_stability/<str:chem_id>/del/<int:id>', views.liver_stability_delete, name='liver_stability_delete'),
    path('chemicals/<str:target>/CYPInhibition/<str:chem_id>/del/<int:id>', views.cyp_inhibition_delete, name='cyp_inhibition_delete'),
# result update
    path('chemicals/<str:target>/pharmacokinetic/<str:chem_id>/update/<int:id>', views.pharmacokinetic_update, name='pharmacokinetic_update'),
    path('chemicals/<str:target>/cytotoxicity/<str:chem_id>/update/<int:id>', views.cytotoxicity_update, name='cytotoxicity_update'),
    path('chemicals/<str:target>/liverstability/<str:chem_id>/update/<int:id>', views.liver_stability_update, name='liver_stability_update'),
    path('chemicals/<str:target>/cyp/<str:chem_id>/update/<int:id>', views.cyp_inhibition_update,name='cyp_inhibition_update'),

    path('chemicals/<str:target>/modeling/<str:chem_id>/', views.schrodinger_model_list, name='schrodinger_model_list'),
    path('chemicals/<str:target>/modeling/<str:chem_id>/add/', views.schrodinger_model_add, name='schrodinger_model_add'),
    path('chemicals/<str:target>/modeling/<str:chem_id>/upload/', views.schrodinger_model_upload, name='schrodinger_model_upload'),
# result invitro add
    path('chemicals/<str:target>/cck/<str:chem_id>/add/',views.cck_add, name='cck_add'),
    path('chemicals/<str:target>/wb/<str:chem_id>/add/',views.wb_add, name='wb_add'),
    path('chemicals/<str:target>/in_target/<str:chem_id>/add/',views.in_target_add, name='in_target_add'),
    path('chemicals/<str:target>/other/<str:chem_id>/add/',views.other_add, name='other_add'),
# result invitro delete
    path('chemicals/<str:target>/cck/<str:chem_id>/del/<int:id>', views.cck_delete, name='cck_delete'),
    path('chemicals/<str:target>/wb/<str:chem_id>/del/<int:id>', views.wb_delete, name='wb_delete'),
    path('chemicals/<str:target>/it/<str:chem_id>/del/<int:id>', views.in_target_delete, name='in_target_delete'),
    path('chemicals/<str:target>/ot/<str:chem_id>/del/<int:id>', views.other_delete, name='other_delete'),
# result invitro update
    path('chemicals/<str:target>/cck/<str:chem_id>/update/<int:id>', views.cck_update, name='cck_update'),
    path('chemicals/<str:target>/wb/<str:chem_id>/update/<int:id>', views.wb_update, name='wb_update'),
    path('get_columns/<str:table_name>/', views.get_columns, name='get_columns'),
    path('save_table_data/', views.save_table_data, name='save_table_data' ),
# FDA
    path('chemicals/<str:target>/FDA/<str:chem_id>', views.FDA_result_view, name='FDA_result_view'),
    path('chemicals/result_fda/<str:target>/', views.upload_fda_result, name='upload_fda_result'),
    path('chemicals/<str:target>/<str:chem_id>', views.download_file, name='download_file'),
    path('delete_file/<int:id>/<str:target>/', views.delete_file, name='delete_file'),
    path('chemicals/<str:target>/FDA_add/<str:chem_id>', views.FDA_result_add, name='FDA_result_add'),
    path('chemicals/<str:target>/FDA_del/<str:chem_id>/<int:id>', views.FDA_delete, name='FDA_delete'),
    path('chemicals/<str:target>/FDA_update/<str:chem_id>/<int:id>', views.FDA_update, name='FDA_update'),
]
if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
