from django.shortcuts import render

# Create your views here.
from django.shortcuts import render, redirect, get_object_or_404
from .models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel, LiverMicrosomalStability, CYPInhibition
from .forms import ChemicalForm, ChemicalUploadForm, PharmacokineticForm, CytotoxicityForm, SchrödingerModelForm, SchrödingerModelUploadForm, LiverMicrosomalStabilityForm, CYPInhibitionForm
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.decorators import login_required
# R2Py
import subprocess
import pandas as pd
# RDKIT
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from io import BytesIO
from django.core.files.base import ContentFile
import logging
import csv
from django.http import JsonResponse
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
    sort_by = request.GET.get('sort_by', 'chem_id')
    sort_order = request.GET.get('sort_order', 'asc')

    if sort_order == 'asc':
        chemicals = Chemical.objects.filter(target=target).order_by(sort_by)
    else:
        chemicals = Chemical.objects.filter(target=target).order_by('-' + sort_by)

    return render(request, 'chemicals/target.html', {
        'chemicals': chemicals,
        'target': target,
        'sort_by': sort_by,
        'sort_order': sort_order,
    })
    # chemicals = Chemical.objects.filter(target__iexact=target)
    # if not chemicals.exists():
    #     logger.debug(f'No chemicals found for target: {target}')
    # else:
    #     logger.debug(f'Found chemicals for target {target}: {chemicals}')
    # return render(request, 'chemicals/target.html', {'chemicals': chemicals, 'target': target,})
@login_required
def delete_selected_chems(request, target):
    if request.method == 'POST':
        selected_chems = request.POST.getlist('selected_chems[]')
        Chemical.objects.filter(pk__in=selected_chems).delete()
        return JsonResponse({'status': 'success'})
    return JsonResponse({'status': 'error'}, status=400)
@login_required
def chemical_new_view(request, target):
    if request.method == 'POST':
        form = ChemicalForm(request.POST, request.FILES)
        if form.is_valid():
            chemical = form.save(commit=False)
            chemical.target = target
            chemical.cLogP = calculate_cLogP(chemical.smiles)
            image_data = generate_image(chemical.smiles)
            if image_data:
                 chemical.image.save(f'{chemical.chem_id}.png', ContentFile(image_data), save=False)
            chemical.save()
            logger.debug(f'Chemical saved: {chemical}')
            return redirect('target_view', target=target)
        else:
            logger.error(f'Form is not valid: {form.errors}')
    else:
        form = ChemicalForm()
    logger.debug(f'New chemical form for {target}')
    return render(request, 'chemicals/chemical_form.html', {'form': form, 'target': target})
def calculate_cLogP(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Descriptors.MolLogP(mol)
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
def chemical_edit_view(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = ChemicalForm(request.POST, instance=chemical)
        if form.is_valid():
            chemical = form.save(commit=False)
            chemical.cLogP = calculate_cLogP(chemical.smiles)
            image_data = generate_image(chemical.smiles)
            if image_data:
               chemical.image.save(f'{chemical.chem_id}.png', ContentFile(image_data), save=False)
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
    liver_stabilities = LiverMicrosomalStability.objects.filter(chemical=chemical)
    cyp_inhibitions = CYPInhibition.objects.filter(chemical=chemical)
    return render(request, 'chemicals/pharmacokinetic_list.html', {
        'chemical': chemical,
        'pharmacokinetics': pharmacokinetics,
        'cytotoxicities': cytotoxicities,
        'liver_stabilities': liver_stabilities,
        'cyp_inhibitions': cyp_inhibitions
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
                    QPlogS=row.get('QPlogS', 0),
                    QPlogHERG=row.get('QPlogHERG', 0),
                    QPPCaco=row.get('QPPCaco', 0),
                    QPlogBB=row.get('QPlogBB', 0),
                    QPPMDCK=row.get('QPPMDCK', 0),
                    Metab=row.get('Metab', 0),
                    QPlogKhsa=row.get('QPlogKhsa', 0),
                    HOralAbs=row.get('HOralAbs', 0),
                    PerHOralAbs=row.get('PerHOralAbs', 0),
                )
                model.save()
            return redirect('schrodinger_model_list', target=target, chem_id=chem_id)
    else:
        form = SchrödingerModelUploadForm()
    return render(request, 'chemicals/schrodinger_model_upload.html', {'form': form, 'chemical': chemical})

@login_required
def liver_stability_add(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = LiverMicrosomalStabilityForm(request.POST)
        if form.is_valid():
            liver_stability = form.save(commit=False)
            liver_stability.chemical = chemical
            liver_stability.save()
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = LiverMicrosomalStabilityForm()
    return render(request, 'chemicals/liver_stability_form.html', {'form': form, 'chemical': chemical})

@login_required
def cyp_inhibition_add(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = CYPInhibitionForm(request.POST)
        if form.is_valid():
            cyp_inhibition = form.save(commit=False)
            cyp_inhibition.chemical = chemical
            cyp_inhibition.save()
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = CYPInhibitionForm()
    return render(request, 'chemicals/cyp_inhibition_form.html', {'form': form, 'chemical': chemical})

def run_r_script_and_get_results(patient_id):
    # R 스크립트 파일 읽기
    with open("/Users/sangjoonshin/Desktop/231218_cosine_similarity.R", "r") as file:
        r_script = file.readlines()

    # selected_patient_id 설정
    patient_id_line = f'selected_patient_id <- "{patient_id}"\n'
    r_script.insert(0, patient_id_line)  # 스크립트 시작 부분에 삽입

    # 수정된 R 스크립트 파일 임시 저장
    with open("/Users/sangjoonshin/Desktop/temp_cosine_jy.R", "w") as file:
        file.writelines(r_script)

    # 수정된 R 스크립트 실행
    subprocess.run(["/usr/local/bin/Rscript", "/Users/sangjoonshin/Desktop/temp_cosine_jy.R"])

    # 결과 파일 읽기
    results_df = pd.read_csv("/Users/sangjoonshin/Desktop/top_10_cell_lines.csv")
    return results_df
def patient_input(request):
    context = {}
    if request.method == 'POST':
        patient_id = request.POST.get('patient_id')
        results_df = run_r_script_and_get_results(patient_id)
        context['results'] = results_df.to_html()  # DataFrame을 HTML로 변환
        context['patient_id'] = patient_id
    return render(request, 'r.html', context)