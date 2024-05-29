from django import forms
from .models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel,LiverMicrosomalStability


class DateInput(forms.DateInput):
    input_type = 'date'


class ChemicalForm(forms.ModelForm):
    class Meta:
        model = Chemical
        fields = ['chem_id', 'smiles', 'image', 'MW', 'cLogP', 'TPSA','H_donors', 'H_acceptors', 'lipinski' ]

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

class SchrödingerModelForm(forms.ModelForm):
    class Meta:
        model = SchrödingerModel
        fields = ['field_1', 'field_2', 'field_3', 'field_4', 'field_5', 'field_6', 'field_7', 'field_8', 'field_9', 'field_10', 'field_11', 'field_12', 'field_13', 'field_14', 'field_15', 'field_16', 'field_17', 'field_18', 'field_19', 'field_20']

class SchrödingerModelUploadForm(forms.Form):
    file = forms.FileField()

class LiverMicrosomalStabilityForm(forms.ModelForm):
    class Meta:
        model = LiverMicrosomalStability
        fields = ['date', 'mouse', 'rat', 'human']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }