from django.urls import path
from .views import create_table_view

urlpatterns = [
    path('data/create-table/', create_table_view, name='create_table_view' ),
]
