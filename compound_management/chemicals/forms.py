from django import forms
from django.forms import ClearableFileInput, formset_factory
from .models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel, LiverMicrosomalStability, CYPInhibition, \
     CCK_assay, invtro_Image, Western_blot, Target_Inhibition, other_asssay, in_vivo
class DateInput(forms.DateInput):
    input_type = 'date'

class ChemicalForm(forms.ModelForm):
    class Meta:
        model = Chemical
        fields = ['chem_id', 'smiles', 'image', 'MW', 'cLogP', 'TPSA' , 'H_donors', 'H_acceptors', 'lipinski','user']

class ChemicalUploadForm(forms.Form):
    file = forms.FileField()

class PharmacokineticForm(forms.ModelForm):
    class Meta:
        model = Pharmacokinetic
        fields = ['date', 'max_concentration', 'tmax', 'AUC', 't_half', 'Vss', 'F', 'CL', 'Route','user']
        widgets = {
            'date': DateInput(attrs={'type': 'date'}),  # DateInput 위젯 사용
        }

class CytotoxicityForm(forms.ModelForm):
    class Meta:
        model = Cytotoxicity
        fields = ['date', 'VERO', 'HFL_1', 'L929', 'NIH_3T3', 'CHO_K1', 'user']
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
        fields = ['date', 'mouse', 'rat', 'human' , 'user']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class CYPInhibitionForm(forms.ModelForm):
    class Meta:
        model = CYPInhibition
        fields = ['date', 'cyp_1a2', 'cyp_2c9', 'cyp_2c19', 'cyp_2d6', 'cyp_3a4' , 'user']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class cckForm(forms.ModelForm):
    class Meta:
        model = CCK_assay
        fields = ['date', 'user',  'cell', 'IC50' , 'Out' , 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }
class wbForm(forms.ModelForm):
    class Meta:
        model = Western_blot
        fields = ['date', 'user' , 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class intargetForm(forms.ModelForm):
    class Meta:
        model = Target_Inhibition
        fields = ['date', 'user', 'vitro_at_10' , 'Taret_IC50' , 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class otherForm(forms.ModelForm):
    class Meta:
        model = other_asssay
        fields = ['date', 'user', 'title' , 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

# class invtroimgForm(forms.ModelForm):
#     images = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True}))
#     class Meta:
#         model = invtro_Image
#         fields = ['image']

class in_vivoForm(forms.ModelForm):
    class Meta:
        model = in_vivo
        fields = ['user', 'start_date', 'end_date', 'cell', 'does', 'solvent', 'inject_date', 'group', 'comment', 'category']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }
class CustomClearableFileInput(ClearableFileInput):
    allow_multiple_selected = True

    def __init__(self, attrs=None):
        super().__init__(attrs)
        if attrs is None:
            attrs = {}
        attrs.update({'multiple': True})

    def value_from_datadict(self, data, files, name):
        if self.allow_multiple_selected:
            return files.getlist(name)
        return files.get(name)

class InvtroimgForm(forms.ModelForm):
    images = forms.ImageField(widget=CustomClearableFileInput)
    class Meta:
        model = invtro_Image
        fields = ['images']

