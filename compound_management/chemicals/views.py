import io
import json
import re

import chardet
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
from django.forms import modelformset_factory, model_to_dict
from django.shortcuts import render, redirect, get_object_or_404
from django.utils.text import capfirst
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_POST
from django.contrib.contenttypes.models import ContentType
from .models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel, LiverMicrosomalStability, CYPInhibition, \
    CCK_assay, invtro_Image, Western_blot, Target_Inhibition, other_asssay, in_vivo, Favorite, FDA_result, Document, \
    result_document
from .forms import ChemicalForm, ChemicalUploadForm, PharmacokineticForm, CytotoxicityForm, SchrödingerModelForm, \
    SchrödingerModelUploadForm, LiverMicrosomalStabilityForm, CYPInhibitionForm, cckForm, wbForm, \
    intargetForm, otherForm, in_vivoForm, FDA_Form, DocumentForm, FDA_UploadForm, InvtroimgForm
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.decorators import login_required
# RDKIT
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from io import BytesIO
from django.core.files.base import ContentFile
import logging
import csv
from django.http import JsonResponse, HttpResponse, HttpResponseNotFound
import os
from users.models import DatabaseList


# 로거 설정
logger = logging.getLogger(__name__)


@login_required
def home_view(request):
    return render(request, 'chemicals/home.html')


@login_required
def target_view(request, target):
    sort_by = request.GET.get('sort_by', 'chem_id')
    sort_order = request.GET.get('sort_order', 'asc')
    favorite_items = Favorite.objects.filter(user=request.user).values_list('item_id',flat=True)

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
        'favorite_items': favorite_items,
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
            if chem.image and chem.image.path and os.path.exists(chem.image.path):  # 이미지가 존재하는지 확인
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

            if smiles != new_value:
                # 파일이 존재하는지 확인하고, 파일 경로가 있을 때만 삭제
                if chemical.image and hasattr(chemical.image, 'path') and os.path.exists(chemical.image.path):
                    os.remove(chemical.image.path)

                # 새로운 이미지를 생성하고 저장
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

        if chemical.image and os.path.isfile(chemical.image.path) : #경로에 있는 파일이 실제파일인지 아닌지
            os.remove(chemical.image.path) # 주어진 경로의 파일 삭제
            # image_folder = os.path.dirname(chemical.image.path) # 같은 경로 반환
            # # png_file_path = os.path.join(image_folder, f"{chemical.chem_id}.png") #같은 경로의 png파일
            # # if os.path.isfile(png_file_path):
            # #     os.remove(png_file_path)

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
            print(reader.fieldnames)
            for row in reader:
                smiles =  row.get('smiles')
                MW = row.get('MW')
                chem_id = row.get('chem_id')
                user = request.user.userID
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


db_names = DatabaseList.objects.all() #쿼리셋
db_list = list(db_names.values_list('name', flat=True)) #리스트

@login_required
def pharmacokinetic_list(request, target, chem_id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)

    cck_assay = CCK_assay.objects.filter(chemical=chemical)
    wb = Western_blot.objects.filter(chemical=chemical)
    in_target = Target_Inhibition.objects.filter(chemical=chemical)
    others = other_asssay.objects.filter(chemical=chemical)

    pharmacokinetics = Pharmacokinetic.objects.filter(chemical=chemical)
    cytotoxicities = Cytotoxicity.objects.filter(chemical=chemical)
    liver_stabilities = LiverMicrosomalStability.objects.filter(chemical=chemical)
    cyp_inhibitions = CYPInhibition.objects.filter(chemical=chemical)

    invivo = in_vivo.objects.filter(chemical=chemical)

    results = {} # 데이터자동화 3
    for db in db_names:
        if db.id >= 10:
            try:
                model = apps.get_model('data', db.name) #모델 가져옴
            except LookupError: # 모델 없을때 예외처리
                print(f"'{db.name}' does not exist.")
                continue  # 오류시, 다음 루프로 건너뜀
            queryset = model.objects.filter(chemical=chemical)
            results[db.id] = {
                'name': db.name,  # db.name 추가
                'data': list(queryset.values())  # 쿼리셋을 리스트로 변환
            }

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
        'invivo': invivo,
        'results': results,
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
    #chemical = get_object_or_404(Chemical, chem_id=chem_id)
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
                invtro_Image.objects.create(CCK_assay=CCK_assay, image=image) # CCK_Image 모델의 새로운 객체를 생성하고 데이터베이스에 저장

            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = cckForm()
    return render(request, 'chemicals/cck_assay_form.html', {
        'form': form,
        'chem_id': chem_id,
        'target': target,
    })


@login_required
def cck_delete(request, target, chem_id, id):
    cck_images = invtro_Image.objects.filter(CCK_assay_id=id)

    # CCK_assay 인스턴스 가져오기
    CCKassay = get_object_or_404(CCK_assay, id=id)

    if request.method == 'POST':
        # 관련 이미지 파일 삭제
        for image in cck_images:
            if image.image and os.path.exists(image.image.path):
                os.remove(image.image.path)

        # CCK_assay 인스턴스 삭제
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
                invtro_Image.objects.create(Western_blot=Western_blot, image=image)
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
    wbimages = invtro_Image.objects.filter(Western_blot=id)
    if request.method == 'POST':
        for images in wbimages:
            if images.image and os.path.exists(images.image.path):
                os.remove(images.image.path)
        Westernblot.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)

@login_required
def wb_update(request, target, chem_id ,id):
    wb = get_object_or_404(Western_blot, id=id)
    if request.method == 'POST':
        form = wbForm(request.POST, instance=wb)
        if form.is_valid() :
            wb = form.save(commit=False)
            wb.save()
            return redirect('pharmacokinetic_list', target=target, chem_id=chem_id)
    else:
        form = wbForm(instance=wb)
    return render(request, 'chemicals/wb_form.html', {'form': form, 'target': target, 'id':id, 'chem_id': chem_id})

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
                invtro_Image.objects.create(Target_Inhibition=Target_Inhibition, image=image)
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
    in_image = invtro_Image.objects.filter(Target_Inhibition=id)
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
                invtro_Image.objects.create(other_asssay=other, image=image)
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
    other_image = invtro_Image.objects.filter(other_asssay=id)
    other = get_object_or_404(other_asssay, id=id)
    if request.method == 'POST' :
        for images in other_image :
            if images.image and os.path.exists(images.image.path) :
                os.remove(images.image.path)
        other.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)

def sanitize_attribute_name(name):
    # 유효한 문자만 허용 (영문자, 숫자, 하이픈, 밑줄)
    return re.sub(r'[^a-zA-Z0-9-_]', '', name)

@login_required
@csrf_exempt
def result_add(request, target):
    chemicals = Chemical.objects.filter(target=target).order_by('-datetime')
    db_names = DatabaseList.objects.all()
    sanitized_chemicals = []
    for chemical in chemicals:
        sanitized_chemical = {
            'chem_id': sanitize_attribute_name(chemical.chem_id),
            'MW': chemical.MW,
            'cLogP': chemical.cLogP,
            'TPSA': chemical.TPSA,
        }
        sanitized_chemicals.append(sanitized_chemical)

    if request.method == 'POST':
        selected_db = request.POST.get('dbOption')
        data_to_save = request.POST.get('data_field')  # 저장할 데이터
        try:
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
        'chemicals': sanitized_chemicals ,
        'db_names': db_names,
    })

@login_required
def get_columns(request, table_name): # 다중 data 업로드 시, 변경data 컬럼을 가져오는 코드
    try:
        model = apps.get_model('chemicals', table_name) # 기본 앱
    except LookupError:
        try:
            model = apps.get_model('data', table_name) # 동적 생성 앱
        except LookupError:
            return JsonResponse({'fields': []})

    fields = []
    for field in model._meta.get_fields():
        if isinstance(field, models.ManyToOneRel): # 관계형 필드 (ManyToOneRel) 패스
            continue
        fields.append({
            'name': field.name,
            'help_text': getattr(field, 'help_text', '') # help_text 속성 없는 예외 처리
        })

    return JsonResponse({'fields': fields})

@login_required
@csrf_exempt
def save_table_data(request):
    if request.method == 'POST':
        try:
            data = request.POST
            files = request.FILES.getlist('images')  # 여러 파일 가져오기

            for key in data:
                if key.startswith('row'):
                    row = json.loads(data[key])  # 각 row는 JSON으로 변환
                    db_name = row.pop('db_name')  # 'db_name' 추출
                    try:
                        model = apps.get_model('chemicals', db_name)  # 기본 앱에서 모델 가져오기
                    except LookupError:
                        try:
                            model = apps.get_model('data', db_name)  # 동적 생성 앱에서 모델 가져오기
                        except LookupError:
                            return JsonResponse({'status': 'error', 'message': f'Model {db_name} not found'}, status=400)

                    if model:
                        chemical_id = row.pop('chemical', None)
                        if chemical_id:
                            try:
                                chemical_instance = Chemical.objects.get(chem_id=chemical_id)
                                row['chemical'] = chemical_instance
                            except Chemical.DoesNotExist:
                                return JsonResponse({'status': 'error', 'message': f'Chemical {chemical_id} not found'}, status=400)

                        # 모델 인스턴스 생성 및 저장
                        instance = model.objects.create(**row)

                        # 동적으로 외래키 필드에 이미지 저장
                        related_name = db_name  # db_name을 사용하여 외래키 필드명 생성
                        if hasattr(invtro_Image, related_name):
                            for image in files:
                                image_instance = invtro_Image(image=image)
                                setattr(image_instance, related_name, instance)  # 외래키 필드를 정확히 설정
                                image_instance.save()
                        else:
                            print(f"Invalid foreign key field: {related_name}")

            return JsonResponse({'status': 'success'}, status=201)
        except Exception as e:
            print("Error:", e)
            return JsonResponse({'status': 'error', 'message': str(e)}, status=400)
    else:
        return JsonResponse({'status': 'invalid request'}, status=400)

@login_required
@require_POST
def toggle_favorite(request, chem_id):
    print(f"Toggle favorite called with chem_id: {chem_id}")  # 이 줄을 추가
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    favorite, created = Favorite.objects.get_or_create(user=request.user, item=chemical)

    if not created:
        favorite.delete()
        is_favorite = False
    else:
        is_favorite = True

    return JsonResponse({'success': True, 'is_favorite': is_favorite})


@login_required
def FDA_result_view(request, target, chem_id, period=None ):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    fda = FDA_result.objects.filter(chemical=chemical)
    # query = FDA_result.objects.filter(target=target)

    if period:
        fda = fda.filter(period=period)

    FDA = fda

    return render(request, 'chemicals/FDA_result.html', {
        'chemical': chemical,
        'target': target,
        'FDA': FDA,
    })

@login_required
def FDA_result_add (request, target, chem_id ):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    if request.method == 'POST':
        form = FDA_Form (request.POST)
        if form.is_valid():
            fda = form.save(commit=False)
            fda.chemical = chemical
            fda.save()
            return redirect('FDA_result_view', target=target, chem_id=chem_id)
    else:
        FDA_Form()
    return render(request, 'chemicals/FDA_result_form.html', {
        'chemical': chemical,
        'target': target,
        'chem_id': chem_id,
    })

@login_required
def FDA_delete(request, target, chem_id ,id):
    fda = get_object_or_404(FDA_result, id=id)
    if request.method == 'POST'  :
        fda.delete()
        return JsonResponse({'success': True, 'id': id})
    return JsonResponse({'success': False, 'error': 'Invalid request'}, status=400)

@login_required
def FDA_update(request, target, chem_id ,id):
    fda = get_object_or_404(FDA_result, id=id)
    if request.method == 'POST':
        form = FDA_Form(request.POST, instance=fda)
        if form.is_valid():
            fda = form.save(commit=False)
            fda.save()
            return redirect('FDA_result_view', target=target, chem_id=chem_id)
    else:
        form = FDA_Form(instance=fda)
    return render(request, 'chemicals/FDA_result_form.html', {'form': form, 'target': target, 'id':id, 'chem_id': chem_id})

@login_required
def download_file(request, chem_id, target): # target 연결
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    documents = chemical.documents.filter(chemical=chemical)

    if request.method == 'POST':
        form = DocumentForm(request.POST, request.FILES)
        if form.is_valid():
            document = form.save(commit=False)
            document.chemical = chemical
            document.save()
            return redirect('download_file',chem_id=chem_id, target=target )
    else:
        form = DocumentForm()

    return render(request, 'chemicals/download_file.html', {
        'documents': documents,
        'target': target,
        'chemical': chemical,
        'form': form
    })

@login_required
def result_file(request, db, id): # result 연결 generic foreignkey 사용
    content_type = ContentType.objects.get(model=db.lower())
    model_class = content_type.model_class()
    product = get_object_or_404(model_class, id=id)

    if request.method == 'POST':
        if 'file' in request.FILES:
            file = request.FILES['file']
            result_document.objects.create(content_object=product, file=file)
            return redirect('result_file', db=db, id=id)
    files = product.files.all()

    return render(request, 'chemicals/result_file.html', {
        'product': product,
        'files': files,
        'db': db,
    })

def delete_file(request,id,target):
    print(target)
    if target == "FDA" :
        documents = get_object_or_404(Document, id=id)
        if request.method == 'POST':
            documents.file.delete()
            documents.delete()
            return JsonResponse({'status': 'success'})
    else :
        file = get_object_or_404(result_document, id=id)
        if request.method == 'POST':
            file.file.delete()
            file.delete()
            return JsonResponse({'status': 'success'})
    return JsonResponse({'status': 'fail'}, status=400)


def upload_fda_result(request,target):
    if request.method == 'POST':
        form = FDA_UploadForm(request.POST, request.FILES)
        if form.is_valid():
            csv_file = request.FILES['file']
            raw_data = csv_file.read()
            # 인코딩 감지
            result = chardet.detect(raw_data)
            encoding = result['encoding']
            print(f"Detected encoding: {encoding}")

            data_set = raw_data.decode(encoding)
            io_string = io.StringIO(data_set)
            reader = csv.DictReader(io_string)
            for row in reader:
                chemical_name=row.get('Product') # 정확한 값을 가져옴
                try:
                    # for row in reader:
                    #     smiles = row['smiles']
                    #     MW = row['MW']
                    #     chem_id = row['chem_id']
                    #     try:
                    #         MW_value = float(MW) if MW else 0
                    #         chemical = Chemical(
                    #             chem_id=chem_id,
                    #             smiles=smiles,
                    #             target=target,
                    #             MW=MW_value,
                    #             cLogP=calculate_cLogP(smiles)
                    #         )
                    #         image_data = generate_image(smiles)
                    #         if image_data:
                    #             chemical.image.save(f'{chemical.chem_id}.png', ContentFile(image_data), save=False)
                    #      chemical.save()

                    chemical = Chemical.objects.get(chem_id=chemical_name)
                    FDA_result.objects.create(
                        chemical=chemical,
                        tmax= row.get('Tmax_1'),
                        max_concentration=row.get('Cmax_1'),
                        AUC=row.get('AUC_1'),
                        t_half=row.get('T1/2_1'),
                        period='1',
                        user=request.user,
                    )
                    FDA_result.objects.create(
                        chemical=chemical,
                        tmax=row.get('Tmax_2'),
                        max_concentration=row.get('Cmax_2'),
                        AUC=row.get('AUC_2'),
                        t_half=row.get('T1/2_2'),
                        period='2',
                        user=request.user
                    )
                except Chemical.DoesNotExist:
                    print(f"Chemical '{chemical_name}' does not exist. Skipping this entry.")
            return redirect('target_view', target=target)
    else:
        form = FDA_UploadForm
    return render(request, 'chemicals/FDA_upload.html', {'form': form, 'target': target})

@login_required # 미완성 : img list
def result_img(request, target, chem_id,id):
    chemical = get_object_or_404(Chemical, chem_id=chem_id)
    image = chemical.invtro_Image.filter(chemical=chemical)

    if request.method == 'POST':
        form = InvtroimgForm(request.POST, request.FILES)
        if form.is_valid():
            image = form.save(commit=False)
            image.chemical = chemical
            image.save()
            return redirect('result_img', chem_id=chem_id, target=target)
    else:
        form = DocumentForm()

    return render(request, 'chemicals/result_img.html', {
        'image': image,
        'target': target,
        'chemical': chemical,
        'form': form
    })
# def ajax_pagination(request):
#     table_id = request.GET.get('table_id')  # 요청된 테이블 ID (또는 유형)
#     page = request.GET.get('page', 1)
#     per_page = 10
#
#     if table_id == 'table1':
#         items = Pharmacokinetic.objects.all()  # 첫 번째 테이블의 데이터
#     elif table_id == 'table2':
#         items = Cytotoxicity.objects.all()  # 두 번째 테이블의 데이터
#     elif table_id == 'table3':
#         items = LiverMicrosomalStability.objects.all()  # 세 번째 테이블의 데이터
#     elif table_id == 'table4':
#         items = CYPInhibition.objects.all()  # 네 번째 테이블의 데이터
#
#     paginator = Paginator(items, per_page)
#     page_obj = paginator.get_page(page)
#
#     data = {
#         'items': list(page_obj.object_list.values()),
#         'has_next': page_obj.has_next(),
#         'has_previous': page_obj.has_previous(),
#         'page_number': page_obj.number,
#         'num_pages': paginator.num_pages,
#     }
#     return JsonResponse(data)