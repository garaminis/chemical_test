from django.contrib import messages
from django.contrib.auth.hashers import check_password
from django.core.paginator import Paginator
from django.forms import modelformset_factory
from django.shortcuts import render

# Create your views here.
from django.shortcuts import render, redirect, get_object_or_404
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_protect
from django.views.decorators.http import require_POST

from .models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel, LiverMicrosomalStability, CYPInhibition, \
    User, CCK_assay, invtro_Image, Western_blot, Target_Inhibition, other_asssay
from .forms import ChemicalForm, ChemicalUploadForm, PharmacokineticForm, CytotoxicityForm, SchrödingerModelForm, \
    SchrödingerModelUploadForm, LiverMicrosomalStabilityForm, CYPInhibitionForm, UserForm, cckForm, wbForm, \
    intargetForm, otherForm
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
import os
from django.conf import settings
from django.http import HttpResponse

# 로거 설정
logger = logging.getLogger(__name__)

def register_view(request):
    if request.method == 'POST':
        form = UserForm(request.POST)
        if form.is_valid():
            form.save()
            return redirect('login')
    else:
        form = UserForm()
    return render(request, 'chemicals/register.html', {'form': form}) # 오류 렌더링

def login_view(request):
    if request.method == 'POST':
        userID = request.POST['userID']     #userID = request.POST.get('userID')도 가능
        password = request.POST['password']
        user = authenticate(request, username=userID, password=password)
        print(user)
        if user is not None:
            login(request, user)
            return redirect('home')
        else:
            messages.error(request, '로그인에 실패하였습니다.')
    return render(request, 'chemicals/login.html')

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

    paginator = Paginator(chemicals, 10)  # 갯수 정해서 보여줌
    page_number = request.GET.get('page') #get요청된 페이지 번호
    page_obj = paginator.get_page(page_number) # 해당 번호에 맞는 페이지 가져옴

    start_index = (page_obj.number - 1) // 10 * 10 + 1
    end_index = min(page_obj.number + 9, paginator.num_pages)
    page_range = range(start_index, end_index + 1)

    return render(request, 'chemicals/target.html', {
        'chemicals': page_obj,
        # 'chemicals': chemicals,
        'target': target,
        'sort_by': sort_by,
        'sort_order': sort_order,
        'page_range': page_range,
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
            user = request.POST.get('user')
            chemical.user = user
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

        if os.path.isfile(chemical.image.path) : #경로에 있는 파일이 실제파일인지 아닌지
            os.remove(chemical.image.path) # 주어진 경로의 파일 삭제
        image_folder = os.path.dirname(chemical.image.path) # 같은 경로 반환
        # png_file_path = os.path.join(image_folder, f"{chemical.chem_id}.png") #같은 경로의 png파일
        # if os.path.isfile(png_file_path):
        #     os.remove(png_file_path)

        chemical.delete()
        logger.debug(f'Chemical deleted: {chemical}')
        return redirect('target_view', target=target)
    return render(request, 'chemicals/chemical_confirm_delete.html', {'chemical': chemical, 'target': target})

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
                try:
                    MW_value = float(MW) if MW else 0
                    chemical = Chemical(
                        chem_id=chem_id,
                        smiles=smiles,
                        target=target,
                        MW=MW_value,
                        cLogP=calculate_cLogP(smiles)
                    )
                    image_data = generate_image(smiles)
                    if image_data:
                        chemical.image.save(f'{chemical.chem_id}.png', ContentFile(image_data), save=False)
                    chemical.save()
                except ValueError as e:
                    print(f"Skipping chemical with chem_id {chem_id} due to error: {e}")
        return redirect('target_view', target=target)
    else:
        form = ChemicalUploadForm()
    return render(request, 'chemicals/upload_chemicals.html', {'form': form, 'target': target})

@login_required
def pharmacokinetic_list(request, target, chem_id, other=None):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)

    pharmacokinetics = Pharmacokinetic.objects.filter(chemical=chemical)
    cytotoxicities = Cytotoxicity.objects.filter(chemical=chemical)
    liver_stabilities = LiverMicrosomalStability.objects.filter(chemical=chemical)
    cyp_inhibitions = CYPInhibition.objects.filter(chemical=chemical)

    cck_assay = CCK_assay.objects.filter(chemical=chemical).all()
    wb = Western_blot.objects.filter(chemical=chemical).all()
    in_target = Target_Inhibition.objects.filter(chemical=chemical).all()
    others = other_asssay.objects.filter(chemical=chemical).all()

    return render(request, 'chemicals/pharmacokinetic_list.html', {
        'target': target,
        'chemical': chemical,
        'pharmacokinetics': pharmacokinetics,
        'cytotoxicities': cytotoxicities,
        'liver_stabilities': liver_stabilities,
        'cyp_inhibitions': cyp_inhibitions,
        'cck_assay':  cck_assay,
        'Western_blot': wb,
        'Target_Inhibition': in_target,
        'other': others
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
    return render(request, 'chemicals/pharmacokinetic_form.html', {'form': form, 'chemical': chemical,  'target': target})

@login_required
def pharmacokinetic_delete (request, target, chem_id ,id):
    Pharm = get_object_or_404(Pharmacokinetic, id=id)
    if request.method == 'POST':
        Pharm.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)

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
def cytotoxicity_delete(request, target, chem_id ,id):
    Cyto = get_object_or_404(Cytotoxicity, id=id)
    if request.method == 'POST':
        Cyto.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)

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
def liver_stability_delete(request, target, chem_id ,id):
    liver = get_object_or_404(LiverMicrosomalStability, id=id)
    if request.method == 'POST':
        liver.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)

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

@login_required
def cyp_inhibition_delete(request, target, chem_id ,id):
    cyp = get_object_or_404(CYPInhibition, id=id)
    if request.method == 'POST':
        cyp.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)
@login_required
def cck_add(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)

    if request.method == 'POST':
        form = cckForm(request.POST)

        if form.is_valid() :
            CCK_assay = form.save(commit=False)
            CCK_assay.chemical = chemical
            CCK_assay.save()

            images = request.FILES.getlist('images') #  파일 업로드를 처리.리스트 반환
            for image in images:
                invtro_Image.objects.create(cck_assay=CCK_assay, image=image) # CCK_Image 모델의 새로운 객체를 생성하고 데이터베이스에 저장

            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = cckForm()
    return render(request, 'chemicals/cck_assay_form.html', {
        'form': form,
        'chemical': chemical,
        'target': target,
    })

@login_required
def cck_delete (request, target, chem_id ,id):
    cck_images = invtro_Image.objects.filter(cck_assay=id)
    CCKassay = get_object_or_404(CCK_assay, id=id) #특정 조건에 맞는 객체를 데이터베이스에서 가져옴,일치하는 객체가 없으면 404 에러를 반환
    if request.method == 'POST':
        for images in cck_images:
            if images.image and os.path.exists(images.image.path):
                os.remove(images.image.path)
        # CCKassay.objects.get(id=id).delete()
        CCKassay.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)


@login_required
def wb_add(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = wbForm(request.POST)
        if form.is_valid() :
            Western_blot = form.save(commit=False)
            Western_blot.chemical = chemical
            Western_blot.save()
            images = request.FILES.getlist('images')
            for image in images:
                invtro_Image.objects.create(wb=Western_blot, image=image)
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = wbForm()
    return render(request, 'chemicals/wb_form.html', {
        'form': form,
        'chemical': chemical,
        'target': target,
    })

@login_required
def wb_delete (request, target, chem_id ,id):
    Westernblot = get_object_or_404(Western_blot, id=id)
    wbimages = invtro_Image.objects.filter(wb=id)
    if request.method == 'POST':
        for images in wbimages:
            if images.image and os.path.exists(images.image.path):
                os.remove(images.image.path)
        Westernblot.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)

@login_required
def in_target_add(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = intargetForm(request.POST)
        if form.is_valid() :
            Target_Inhibition = form.save(commit=False)
            Target_Inhibition.chemical = chemical
            Target_Inhibition.save()
            images = request.FILES.getlist('images')
            for image in images:
                invtro_Image.objects.create(in_target=Target_Inhibition, image=image)
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = intargetForm()
    return render(request, 'chemicals/in_target_form.html', {
        'form': form,
        'chemical': chemical,
        'target': target,
    })
@login_required
def in_target_delete(request, target, chem_id ,id):
    in_image = invtro_Image.objects.filter(in_target=id)
    in_target = get_object_or_404(Target_Inhibition, id=id)
    if request.method == 'POST' :
        for images in in_image :
            if images.image and os.path.exists(images.image.path) :
                os.remove(images.image.path)
        in_target.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)

@login_required
def other_add(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = otherForm(request.POST)
        if form.is_valid() :
            other = form.save(commit=False)
            other.chemical = chemical
            other.save()
            images = request.FILES.getlist('images')
            for image in images:
                invtro_Image.objects.create(others=other, image=image)
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = intargetForm()
    return render(request, 'chemicals/other_form.html', {
        'form': form,
        'chemical': chemical,
        'target': target,
    })

@login_required
def other_delete(request, target, chem_id ,id):
    other_image = invtro_Image.objects.filter(others=id)
    other = get_object_or_404(other_asssay, id=id)
    if request.method == 'POST' :
        for images in other_image :
            if images.image and os.path.exists(images.image.path) :
                os.remove(images.image.path)
        other.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)

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

def run_r_SL(SLselected_gene):
    # R 스크립트 파일 읽기
    with open("/Users/sangjoonshin/Downloads/SL/Synthetic_Lethal.R", "r") as file:
        r_script = file.readlines()

    # selected_patient_id 설정
    file_path = f"/Users/sangjoonshin/Downloads/SL/DEGene_{SLselected_gene}.csv"
    SLselected_gene = f'SLselected_gene <- "{SLselected_gene}"\n'
    r_script.insert(0, SLselected_gene)  # 스크립트 시작 부분에 삽입

    # 수정된 R 스크립트 파일 임시 저장
    with open("/Users/sangjoonshin/Downloads/SL/temp_SL.R", "w") as file:
        file.writelines(r_script)

    # 수정된 R 스크립트 실행
    subprocess.run(["/usr/local/bin/Rscript", "/Users/sangjoonshin/Downloads/SL/temp_SL.R"])

    # 결과 파일 읽기

    results_df = pd.read_csv(file_path)
    return results_df
def SLselected_gene_input(request):
    context = {}
    if request.method == 'POST':
        SLselected_gene = request.POST.get('SLselected_gene')
        results_df = run_r_SL(SLselected_gene)
        context['results'] = results_df.to_html()  # DataFrame을 HTML로 변환
        context['SLselected_gene'] = SLselected_gene

        # 두 번째 기능: 이미지 경로 설정
        img_relative_path = f'images/CRISPR_Gene_Effect_{SLselected_gene}.png'
        full_img_path = os.path.join(settings.MEDIA_ROOT, img_relative_path)
        if os.path.exists(full_img_path):
            context['img_path'] = os.path.join(settings.MEDIA_URL, img_relative_path)
        else:
            context['img_path'] = None
    return render(request, 'sl.html', context)

