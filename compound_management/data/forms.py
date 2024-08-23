from django import forms
from django.forms import formset_factory


class ColumnForm(forms.Form):
    column_name = forms.CharField(label='Column Name', max_length=100)
    column_type = forms.ChoiceField(
        label='Column Type',
        choices=[('CharField', 'CharField'), ('IntegerField', 'IntegerField')]
    )
ColumnFormSet = formset_factory(ColumnForm, extra=1)

class TableForm(forms.Form):
    table_name = forms.CharField(label='Table Name', max_length=100)