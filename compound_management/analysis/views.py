import os
from django.conf import settings
# R2Py
import subprocess
import pandas as pd
from django.shortcuts import render

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
    return render(request, 'analysis/r.html', context)

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
    return render(request, 'analysis/sl.html', context)

