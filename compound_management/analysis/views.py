import os

from celery.result import AsyncResult
from django.conf import settings
# R2Py
import subprocess
import pandas as pd
from django.shortcuts import render
# Celery
from analysis.task import run_r_script_task

def run_r_script_and_get_results(patient_id):
    # R 스크립트 파일 읽기
    with open("C:/Users/hasn0737/Desktop/r_study/SL/r_file/231218_cosine_similarity.R", "r", encoding='utf-8') as file:   # C:\Users\hasn0737\Desktop\r_study\easy_r\script3.R
        r_script = file.readlines() # 파일의 모든 데이터를 읽어서 각 줄을 리스트 형태로 반환

    # selected_patient_id 설정
    patient_id_line = f'selected_patient_id <- "{patient_id}"\n'
    r_script.insert(0, patient_id_line)  # 스크립트 시작 부분에 삽입

    # 수정된 R 스크립트 파일 임시 저장
    with open("C:/Users/hasn0737/Desktop/r_study/SL/r_file/temp_cosine_jy.R", "w") as file:
        file.writelines(r_script)

    # 수정된 R 스크립트 실행
    subprocess.run(["C:/Program Files/R/R-4.4.1/bin/Rscript", "C:/Users/hasn0737/Desktop/r_study/SL/r_file/temp_cosine_jy.R"])

    # 결과 파일 읽기
    results_df = pd.read_csv("C:/Users/hasn0737/Desktop/r_study/SL/csv_file/top_10_cell_lines.csv")
    return results_df

def patient_input(request):
    context = {}
    if request.method == 'POST':
        patient_id = request.POST.get('patient_id')
        print(patient_id)
        results_df = run_r_script_and_get_results(patient_id) #함수 결과 반환
        print(results_df)
        context['results'] = results_df.to_html()  # DataFrame을 HTML로 변환
        context['patient_id'] = patient_id
    return render(request, 'analysis/r.html', context)

def run_r_SL(SLselected_gene):
    # R 스크립트 파일 읽기
    with open("C:/Users/hasn0737/Desktop/r_study/SL/SL_gr.R", "r", encoding='utf-8') as file:
        r_script = file.readlines()

    # selected_patient_id 설정
    SLselected_gene = f'SLselected_gene <- "{SLselected_gene}"\n'
    r_script.insert(0, SLselected_gene)  # 스크립트 시작 부분에 삽입

    with open("C:/Users/hasn0737/Desktop/r_study/SL/temp_SL.R", "w") as file:
        file.writelines(r_script)

    #file_path = f"C:/Users/hasn0737/Desktop/r_study/SL/csv_file/{SLselected_gene}.csv"
    # 수정된 R 스크립트 실행
    subprocess.run(
        ["C:/Program Files/R/R-4.4.1/bin/Rscript", "C:/Users/hasn0737/Desktop/r_study/SL/temp_SL.R"],
        encoding="utf-8"
    )

    # 결과 파일 읽기
    #results_df = pd.read_csv(file_path)
    return SLselected_gene
    #return results_df

def SLselected_gene_input(request):
    #run_r_script_task.delay()

    context = {}
    if request.method == 'POST':
        SLselected_gene = request.POST.get('SLselected_gene')
        results_df = run_r_SL(SLselected_gene)
        print (SLselected_gene)
        #context['results'] = results_df.to_html()  # DataFrame을 HTML로 변환
        context['SLselected_gene'] = SLselected_gene

        # 두 번째 기능: 이미지 경로 설정
        img_relative_path = f'analysis/{SLselected_gene}.png'
        full_img_path = os.path.join(settings.MEDIA_ROOT, img_relative_path)
        print(full_img_path)
        if os.path.exists(full_img_path):
            context['img_path'] = os.path.join(settings.MEDIA_URL, img_relative_path)
        else:
            context['img_path'] = None
    return render(request, 'analysis/sl.html', context)

