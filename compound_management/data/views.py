import os

from django.apps import apps
from django.contrib.admin.views.decorators import staff_member_required
from django.core.management import call_command
from django.shortcuts import redirect, render
from django.utils.text import capfirst
from .forms import TableForm, ColumnFormSet
from users.models import DatabaseList
from django.db import models

def create_dynamic_model(table_name, fields):
    attrs = {
        '__module__': 'data.models',
        'Meta': type('Meta', (), {'db_table': table_name, 'app_label': 'data'}),
    }
    for field_name, field_type in fields.items():
        if field_type == 'CharField':
            attrs[field_name] = models.CharField(max_length=255)
        elif field_type == 'IntegerField':
            attrs[field_name] = models.IntegerField()
    # 모델 클래스 생성
    model = type(capfirst(table_name), (models.Model,), attrs)

    # 모델을 전역 네임스페이스에 추가
    apps.get_app_config('data').models[model.__name__.lower()] = model


    # DatabaseList에 새로운 DB 이름 추가
    if not DatabaseList.objects.filter(name=table_name).exists():
        DatabaseList.objects.create(name=table_name)

    models_path = os.path.join('data', 'models.py') # a:추가모드
    with open(models_path, 'a') as models_file:
        models_file.write(f"\n\nclass {model.__name__}(models.Model):\n")
        for field_name, field in fields.items():
            if field == 'CharField':
                models_file.write(f"    {field_name} = models.CharField(max_length=255)\n")
            elif field == 'IntegerField':
                models_file.write(f"    {field_name} = models.IntegerField()\n")
        models_file.write(f"    class Meta:\n")
        models_file.write(f"        db_table = '{table_name}'\n")
        models_file.write(f"        app_label = 'data'\n")

    # admin.py에 모델 등록 추가
    admin_path = os.path.join('data', 'admin.py')
    with open(admin_path, 'a') as admin_file:
        admin_file.write(f"\nadmin.site.register({model.__name__})\n")

    return model

@staff_member_required
def create_table_view(request):
    if request.method == 'POST':
        table_form = TableForm(request.POST)
        column_formset = ColumnFormSet(request.POST, prefix='columns')

        if table_form.is_valid() and column_formset.is_valid():
            table_name = table_form.cleaned_data['table_name'].replace(' ', '_').lower()
            columns = column_formset.cleaned_data

            # 필드 정보 구성
            fields = {column['column_name']: column['column_type'] for column in columns if column['column_name'] and column['column_type']}

            # 동적으로 모델을 생성
            create_dynamic_model(table_name, fields)

            call_command('makemigrations', 'data')
            call_command('migrate', 'data')

        return redirect('home')
    else:
        table_form = TableForm()
        column_formset = ColumnFormSet(prefix='columns')

    return render(request, 'admin/create_table.html', {
        'table_form': table_form,
        'column_formset': column_formset,
    })

