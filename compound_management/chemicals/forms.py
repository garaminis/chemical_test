from django import forms
from .models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel,LiverMicrosomalStability, CYPInhibition


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
        fields = ['QPlogS', 'QPlogHERG', 'QPPCaco', 'QPlogBB', 'QPPMDCK', 'Metab', 'QPlogKhsa', 'HOralAbs', 'PerHOralAbs', 'field_11']

class SchrödingerModelUploadForm(forms.Form):
    file = forms.FileField()

class LiverMicrosomalStabilityForm(forms.ModelForm):
    class Meta:
        model = LiverMicrosomalStability
        fields = ['date', 'mouse', 'rat', 'human']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class CYPInhibitionForm(forms.ModelForm):
    class Meta:
        model = CYPInhibition
        fields = ['date', 'cyp_1a2', 'cyp_2c9', 'cyp_2c19', 'cyp_2d6', 'cyp_3a4']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }