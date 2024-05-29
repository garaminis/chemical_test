from django.shortcuts import render

# Create your views here.
from django.shortcuts import render, redirect, get_object_or_404
from .models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel
from .forms import ChemicalForm, ChemicalUploadForm, PharmacokineticForm, CytotoxicityForm, SchrödingerModelForm, SchrödingerModelUploadForm
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.decorators import login_required

from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from io import BytesIO
from django.core.files.base import ContentFile
import logging
import csv
from django.http import HttpResponse
# 로거 설정
logger = logging.getLogger(__name__)


def register_view(request):
    if request.method == 'POST':
        form = UserCreationForm(request.POST)
        if form.is_valid():
            form.save()
            return redirect('login')
    else:
        form = UserCreationForm()
    return render(request, 'chemicals/register.html', {'form': form})

def logout_view(request):
    logout(request)
    return redirect('login')

@login_required
def home_view(request):
    return render(request, 'chemicals/home.html')

@login_required
def target_view(request, target):
    chemicals = Chemical.objects.filter(target__iexact=target)
    if not chemicals.exists():
        logger.debug(f'No chemicals found for target: {target}')
    else:
        logger.debug(f'Found chemicals for target {target}: {chemicals}')
    return render(request, 'chemicals/target.html', {'chemicals': chemicals, 'target': target})
def calculate_cLogP(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return round(Descriptors.MolLogP(mol), 1)
    return None

def generate_image(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        return buffer.getvalue()
    return None

@login_required
def chemical_new_view(request, target):
    if request.method == 'POST':
        form = ChemicalForm(request.POST, request.FILES)
        if form.is_valid():
            chemical = form.save(commit=False)
            # chemical.target = target
            # chemical.cLogP = calculate_cLogP(chemical.smiles)
            # image_data = generate_image(chemical.smiles)
            # if image_data:
            #     chemical.image.save(f'{chemical.chem_id}.png', ContentFile(image_data), save=False)
            chemical.save()
            logger.debug(f'Chemical saved: {chemical}')
            return redirect('target_view', target=target)
        else:
            logger.error(f'Form is not valid: {form.errors}')
    else:
        form = ChemicalForm()
    logger.debug(f'New chemical form for {target}')
    return render(request, 'chemicals/chemical_form.html', {'form': form, 'target': target})

@login_required
def chemical_edit_view(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = ChemicalForm(request.POST, instance=chemical)
        if form.is_valid():
            chemical = form.save(commit=False)
            # chemical.cLogP = calculate_cLogP(chemical.smiles)
            # image_data = generate_image(chemical.smiles)
            # if image_data:
            #     chemical.image.save(f'{chemical.chem_id}.png', ContentFile(image_data), save=False)
            chemical.save()
            logger.debug(f'Chemical updated: {chemical}')
            return redirect('target_view', target=target)
        else:
            logger.error(f'Form is not valid: {form.errors}')
    else:
        form = ChemicalForm(instance=chemical)
    logger.debug(f'Edit chemical form for {chemical}')
    return render(request, 'chemicals/chemical_form.html', {'form': form, 'target': target})
@login_required
def chemical_delete_view(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        chemical.delete()
        logger.debug(f'Chemical deleted: {chemical}')
        return redirect('target_view', target=target)
    return render(request, 'chemicals/chemical_confirm_delete.html', {'chemical': chemical, 'target': target})
def login_view(request):
    if request.method == 'POST':
        username = request.POST['username']
        password = request.POST['password']
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            return redirect('home')
    return render(request, 'chemicals/login.html')
@login_required
def upload_chemicals(request, target):
    if request.method == 'POST':
        form = ChemicalUploadForm(request.POST, request.FILES)
        if form.is_valid():
            file = request.FILES['file']
            decoded_file = file.read().decode('utf-8').splitlines()
            reader = csv.DictReader(decoded_file)
            for row in reader:
                smiles = row['smiles']
                MW = row['MW']
                chem_id = row['chem_id']
                chemical = Chemical(
                    chem_id=chem_id,
                    smiles=smiles,
                    target=target,
                    MW=float(MW),
                    cLogP=calculate_cLogP(smiles)
                )
                image_data = generate_image(smiles)
                if image_data:
                    chemical.image.save(f'{chemical.chem_id}.png', ContentFile(image_data), save=False)
                chemical.save()
        return redirect('target_view', target=target)
    else:
        form = ChemicalUploadForm()
    return render(request, 'chemicals/upload_chemicals.html', {'form': form, 'target': target})

@login_required
def pharmacokinetic_list(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    pharmacokinetics = Pharmacokinetic.objects.filter(chemical=chemical)
    cytotoxicities = Cytotoxicity.objects.filter(chemical=chemical)
    return render(request, 'chemicals/pharmacokinetic_list.html', {
        'chemical': chemical,
        'pharmacokinetics': pharmacokinetics,
        'cytotoxicities': cytotoxicities
    })

@login_required
def pharmacokinetic_add(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = PharmacokineticForm(request.POST)
        if form.is_valid():
            pharmacokinetic = form.save(commit=False)
            pharmacokinetic.chemical = chemical
            pharmacokinetic.save()
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = PharmacokineticForm()
    return render(request, 'chemicals/pharmacokinetic_form.html', {'form': form, 'chemical': chemical})

@login_required
def cytotoxicity_add(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = CytotoxicityForm(request.POST)
        if form.is_valid():
            cytotoxicity = form.save(commit=False)
            cytotoxicity.chemical = chemical
            cytotoxicity.save()
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = CytotoxicityForm()
    return render(request, 'chemicals/cytotoxicity_form.html', {'form': form, 'chemical': chemical})


@login_required
def schrodinger_model_list(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    models = SchrödingerModel.objects.filter(chemical=chemical)
    return render(request, 'chemicals/schrodinger_model_list.html', {'chemical': chemical, 'models': models})

@login_required
def schrodinger_model_add(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = SchrödingerModelForm(request.POST)
        if form.is_valid():
            model = form.save(commit=False)
            model.chemical = chemical
            model.save()
            return redirect('schrodinger_model_list', target=target, chem_id=chem_id)
    else:
        form = SchrödingerModelForm()
    return render(request, 'chemicals/schrodinger_model_form.html', {'form': form, 'chemical': chemical})

@login_required
def schrodinger_model_upload(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = SchrödingerModelUploadForm(request.POST, request.FILES)
        if form.is_valid():
            csv_file = request.FILES['file']
            decoded_file = csv_file.read().decode('utf-8').splitlines()
            reader = csv.DictReader(decoded_file)
            for row in reader:
                model = SchrödingerModel(
                    chemical=chemical,
                    field_1=row.get('field_1', 0),
                    field_2=row.get('field_2', 0),
                    field_3=row.get('field_3', 0),
                    field_4=row.get('field_4', 0),
                    field_5=row.get('field_5', 0),
                    field_6=row.get('field_6', 0),
                    field_7=row.get('field_7', 0),
                    field_8=row.get('field_8', 0),
                    field_9=row.get('field_9', 0),
                    field_10=row.get('field_10', 0),
                    field_11=row.get('field_11', 0),
                    field_12=row.get('field_12', 0),
                    field_13=row.get('field_13', 0),
                    field_14=row.get('field_14', 0),
                    field_15=row.get('field_15', 0),
                    field_16=row.get('field_16', 0),
                    field_17=row.get('field_17', 0),
                    field_18=row.get('field_18', 0),
                    field_19=row.get('field_19', 0),
                    field_20=row.get('field_20', 0)
                )
                model.save()
            return redirect('schrodinger_model_list', target=target, chem_id=chem_id)
    else:
        form = SchrödingerModelUploadForm()
    return render(request, 'chemicals/schrodinger_model_upload.html', {'form': form, 'chemical': chemical})