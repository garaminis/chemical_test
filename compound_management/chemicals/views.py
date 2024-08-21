import re

from django.apps import apps
from django.contrib import messages, admin
from django.contrib.admin.views.decorators import staff_member_required
from django.core.management import call_command
from django.db import connection, models, connections
from django.db.models.base import ModelBase
from django.urls import path
from django.contrib.auth.hashers import check_password
from django.core.paginator import Paginator
from django.db.models import Case, When, IntegerField, CharField, Value
from django.db.models.functions import Cast, Substr, Length
from django.forms import modelformset_factory
from django.shortcuts import render, redirect, get_object_or_404
from django.utils.text import capfirst
from django.views.decorators.csrf import csrf_exempt

from .models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel, LiverMicrosomalStability, CYPInhibition, \
    CCK_assay, invtro_Image, Western_blot, Target_Inhibition, other_asssay, in_vivo
from .forms import ChemicalForm, ChemicalUploadForm, PharmacokineticForm, CytotoxicityForm, SchrödingerModelForm, \
    SchrödingerModelUploadForm, LiverMicrosomalStabilityForm, CYPInhibitionForm, cckForm, wbForm, \
    intargetForm, otherForm, in_vivoForm, TableForm, ColumnFormSet
from django.contrib.auth import authenticate, login, logout
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
from django.http import JsonResponse, HttpResponse
import os
from django.conf import settings
from users.models import DatabaseList


# 로거 설정
logger = logging.getLogger(__name__)

@login_required
def home_view(request):
    return render(request, 'chemicals/home.html')

# db_names = DatabaseList.objects.all()
# db_names = DatabaseList.objects.values_list('name','group')
@login_required
@csrf_exempt
def result_add(request, target):
    db_names = DatabaseList.objects.all()
    chemicals = Chemical.objects.filter(target=target).order_by('-datetime')

    if request.method == 'POST':
        selected_db = request.POST.get('dbOption')
        print(selected_db)
        data_to_save = request.POST.get('data_field')  # 저장할 데이터
        try:
            # 선택된 데이터베이스 테이블에 해당하는 모델 클래스를 동적으로 가져오기
            DynamicModel = apps.get_model('chemicals', selected_db.capitalize())
            if DynamicModel:
                new_entry = DynamicModel(name=data_to_save)
                new_entry.save()
                print(f"Data saved to {selected_db} table")
            else:
                print('해당 모델을 찾을 수 없습니다.')

        except LookupError:
            print(f"The model for {selected_db} was not found.")
        except Exception as e:
            print(f"An error occurred: {str(e)}")

    return render(request, 'chemicals/result_form.html', {
        'target': target,
        'chemicals': chemicals,
        'db_names': db_names,

    })

@login_required
def target_view(request, target):
    sort_by = request.GET.get('sort_by', 'chem_id')
    sort_order = request.GET.get('sort_order', 'asc')

    def sort_key(value):
        str_value = str(value)
        match = re.match(r'([a-zA-Z]*)([\d.]+)', str_value) # 문자와 숫자를 분리
        if match:
            alpha_part = match.group(1)  # 문자 부분
            numeric_part = match.group(2)  # 숫자 부분
            # 숫자 부분을 float으로 변환
            try:
                numeric_part = float(numeric_part)
            except ValueError:
                numeric_part = 0
            return alpha_part, numeric_part
        return str_value, 0

    if sort_order == 'asc':
        chemicals = sorted(
            Chemical.objects.filter(target=target),
            key=lambda obj: sort_key(getattr(obj, sort_by))
        )
    else:
        chemicals = sorted(
            Chemical.objects.filter(target=target),
            key=lambda obj: sort_key(getattr(obj, sort_by)),
            reverse=True
        )

    # if sort_order == 'asc':
    #     chemicals = Chemical.objects.filter(target=target).order_by(sort_by)
    # else:
    #     chemicals = Chemical.objects.filter(target=target).order_by('-' + sort_by)

    # paginator = Paginator(chemicals, 10)  # 갯수 정해서 보여줌
    # page_number = request.GET.get('page') #get요청된 페이지 번호
    # page_obj = paginator.get_page(page_number) # 해당 번호에 맞는 페이지 가져옴
    #
    # start_index = (page_obj.number - 1) // 10 * 10 + 1
    # end_index = min(page_obj.number + 9, paginator.num_pages)
    # page_range = range(start_index, end_index + 1)

    return render(request, 'chemicals/target.html', {
        # 'chemicals': page_obj,
        # 'chemicals': chemicals,
        'chemicals': chemicals,
        'target': target,
        'sort_by': sort_by,
        'sort_order': sort_order,
        # 'page_range': page_range,
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
        chems_img = Chemical.objects.filter(pk__in=selected_chems)
        for chem in chems_img:
            if os.path.exists(chem.image.path):
                os.remove(chem.image.path)
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

            new_value = form.cleaned_data['smiles']
            chemical_instance = Chemical.objects.get(chem_id=chem_id)
            smiles = chemical_instance.smiles
            if smiles != new_value :
                if os.path.exists(chemical.image.path):
                    os.remove(chemical.image.path)
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

    invivo = in_vivo.objects.filter(chemical=chemical).all()

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
        'other': others,
        'invivo' : invivo,
    })

@login_required
def pharmacokinetic_add(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = PharmacokineticForm(request.POST)
        if form.is_valid():
            pharmacokinetic = form.save(commit=False)
            pharmacokinetic.chemical = chemical
            # Pharmacokinetic 모델이 Chemical 모델과 외래 키(ForeignKey)로 연결되어 있고, 해당 필드가 null을 허용하지 않는다면, 이 설정이 필수적임
            pharmacokinetic.save()
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = PharmacokineticForm()
    return render(request, 'chemicals/pharmacokinetic_form.html', {'form': form, 'target': target,'chem_id': chem_id})
@login_required
def pharmacokinetic_delete (request, target, chem_id ,id):
    Pharm = get_object_or_404(Pharmacokinetic, id=id)
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        Pharm.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)
@login_required
def pharmacokinetic_update(request, target, chem_id ,id):
    Pharm = get_object_or_404(Pharmacokinetic, id=id)
    if request.method == 'POST':
        form = PharmacokineticForm(request.POST, instance=Pharm)
        if form.is_valid():
            Pharm = form.save(commit=False)
            Pharm.save()
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = PharmacokineticForm(instance=Pharm)
    return render(request, 'chemicals/pharmacokinetic_form.html', {'form': form, 'target': target, 'id':id,'chem_id': chem_id})
 # form으로 인자전잘 꼭 해줘야함.

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
    return render(request, 'chemicals/cytotoxicity_form.html', {'form': form, 'target': target,'chem_id': chem_id})

@login_required
def cytotoxicity_delete(request, target, chem_id ,id):
    Cyto = get_object_or_404(Cytotoxicity, id=id)
    if request.method == 'POST':
        Cyto.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)

@login_required
def cytotoxicity_update(request, target, chem_id ,id):
    Cyto = get_object_or_404(Cytotoxicity, id=id)
    if request.method == 'POST':
        form = CytotoxicityForm(request.POST, instance=Cyto)
        if form.is_valid():
            Cyto = form.save(commit=False)
            Cyto.save()
            logger.debug(f'Chemical updated: {Cyto}')
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = CytotoxicityForm(instance=Cyto)
    return render(request, 'chemicals/Cytotoxicity_form.html', {'form': form, 'target': target, 'id':id,'chem_id': chem_id})

@login_required
def invivo_add(request, target, chem_id, category):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = in_vivoForm(request.POST)
        if form.is_valid():
            in_vivo = form.save(commit=False)
            in_vivo.chemical = chemical
            in_vivo.save()
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = CytotoxicityForm()
    return render(request, 'chemicals/in_vivo_form.html', {
        'form': form,
        'target': target,
        'chem_id': chem_id,
        'category': category,
    })

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
    return render(request, 'chemicals/liver_stability_form.html', {'form': form, 'target': target,'chem_id': chem_id})

@login_required
def liver_stability_delete(request, target, chem_id ,id):
    liver = get_object_or_404(LiverMicrosomalStability, id=id)
    if request.method == 'POST':
        liver.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)

@login_required
def liver_stability_update(request, target, chem_id ,id):
    liver = get_object_or_404(LiverMicrosomalStability, id=id)
    if request.method == 'POST':
        form = LiverMicrosomalStabilityForm(request.POST, instance=liver)
        if form.is_valid():
            liver = form.save(commit=False)
            liver.save()
            logger.debug(f'Chemical updated: {liver}')
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = LiverMicrosomalStabilityForm(instance=liver)
    return render(request, 'chemicals/liver_stability_form.html', {'form': form, 'target': target, 'id':id, 'chem_id': chem_id })

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
    return render(request, 'chemicals/cyp_inhibition_form.html', {'form': form, 'target': target,'chem_id': chem_id})

@login_required
def cyp_inhibition_delete(request, target, chem_id ,id):
    cyp = get_object_or_404(CYPInhibition, id=id)
    if request.method == 'POST':
        cyp.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)

@login_required
def cyp_inhibition_update(request, target, chem_id ,id):
    CYP = get_object_or_404(CYPInhibition, id=id)
    if request.method == 'POST':
        form = CYPInhibitionForm(request.POST, instance=CYP)
        if form.is_valid():
            CYP = form.save(commit=False)
            CYP.save()
            logger.debug(f'Chemical updated: {CYP}')
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = CYPInhibitionForm(instance=CYP)
    return render(request, 'chemicals/cyp_inhibition_form.html', {'form': form, 'target': target, 'id':id, 'chem_id': chem_id})
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
        'chem_id': chem_id,
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
def cck_update(request, target, chem_id ,id):
    cck = get_object_or_404(CCK_assay, id=id)
    # cck_images = invtro_Image.objects.filter(cck_assay=id)
    if request.method == 'POST':
        form = cckForm(request.POST, instance=cck)
        #formset = invtro_Image(request.POST, request.FILES, instance=cck)
        if form.is_valid() :
            cck = form.save(commit=False)
            cck.save()
            #formset.save()
            logger.debug(f'Chemical updated: {cck}')
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = cckForm(instance=cck)
    return render(request, 'chemicals/cck_assay_form.html', {'form': form, 'target': target, 'id':id, 'chem_id': chem_id})

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
        'chem_id': chem_id,
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
        'chem_id': chem_id,
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
        'chem_id': chem_id,
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


def create_dynamic_model(table_name, fields):
    attrs = {
        '__module__': 'chemicals.models',
        'Meta': type('Meta', (), {'db_table': table_name, 'app_label': 'chemicals'}),
    }
    for field_name, field_type in fields.items():
        if field_type == 'CharField':
            attrs[field_name] = models.CharField(max_length=255)
        elif field_type == 'IntegerField':
            attrs[field_name] = models.IntegerField()

    # 모델 클래스 생성
    model = type(capfirst(table_name), (models.Model,), attrs)

    # 모델을 전역 네임스페이스에 추가
    globals()[model.__name__] = model

    # DatabaseList에 새로운 DB 이름 추가
    if not DatabaseList.objects.filter(name=table_name).exists():
        DatabaseList.objects.create(name=table_name)

    models_path = os.path.join('chemicals', 'models.py') # a:추가모드
    with open(models_path, 'a') as models_file:
        models_file.write(f"\n\nclass {model.__name__}(models.Model):\n")
        for field_name, field in fields.items():
            if field == 'CharField':
                models_file.write(f"    {field_name} = models.CharField(max_length=255)\n")
            elif field == 'IntegerField':
                models_file.write(f"    {field_name} = models.IntegerField()\n")
        models_file.write(f"    class Meta:\n")
        models_file.write(f"        db_table = '{table_name}'\n")
        models_file.write(f"        app_label = 'chemicals'\n")

    # admin.py에 모델 등록 추가
    admin_path = os.path.join('chemicals', 'admin.py')
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

            call_command('makemigrations', 'chemicals')  # 'chemicals' 앱에 대해 마이그레이션 생성
            call_command('migrate', 'chemicals')

        return redirect('home')
    else:
        table_form = TableForm()
        column_formset = ColumnFormSet(prefix='columns')

    return render(request, 'admin/create_table.html', {
        'table_form': table_form,
        'column_formset': column_formset,
    })
