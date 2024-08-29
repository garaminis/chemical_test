from django.urls import path
from .views import create_table_view
from . import views

urlpatterns = [
    path('data/create-table/', create_table_view, name='create_table_view' ),
]
