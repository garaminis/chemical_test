

# Register your models here.
# admin.py
from django.contrib import admin
from .models import Chemical, Result

admin.site.register(Chemical)
admin.site.register(Result)