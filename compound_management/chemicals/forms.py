from django import forms
from .models import Chemical, Pharmacokinetic, Cytotoxicity

class DateInput(forms.DateInput):
    input_type = 'date'


class ChemicalForm(forms.ModelForm):
    class Meta:
        model = Chemical
        fields = ['chem_id', 'smiles', 'image', 'MW', 'cLogP']

class ChemicalUploadForm(forms.Form):
    file = forms.FileField()

class PharmacokineticForm(forms.ModelForm):
    class Meta:
        model = Pharmacokinetic
        fields = ['date', 'cmax', 'tmax', 'AUC', 't_half', 'Vss', 'Vd', 'BA']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class CytotoxicityForm(forms.ModelForm):
    class Meta:
        model = Cytotoxicity
        fields = ['date', 'VERO', 'HFL_1', 'L929', 'NIH_3T3', 'CHO_K1']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }