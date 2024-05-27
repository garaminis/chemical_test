from django import forms
from .models import Chemical

class ChemicalForm(forms.ModelForm):
    class Meta:
        model = Chemical
        fields = ['chem_id', 'smiles', 'image', 'MW', 'cLogP']

class ChemicalUploadForm(forms.Form):
    file = forms.FileField()