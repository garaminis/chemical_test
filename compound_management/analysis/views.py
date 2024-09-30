import os
import time
from celery.result import AsyncResult
from django.conf import settings
# R2Py
import subprocess
import pyRserve
import pandas as pd
from django.shortcuts import render
# Celery
from analysis.task import run_r_script_task
from django.http import JsonResponse
from django.core.cache import cache


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

    # result = cache.get('r_script_result')
    # print(result)
    # if result is None:
    #     return JsonResponse({'status': 'error', 'message': 'No result found in cache.'})

    # R 스크립트 파일 읽기
    with open("C:/Users/hasn0737/Desktop/r_study/SL/Sl_gr.R", "r", encoding='utf-8') as file:
        r_script = file.readlines()

    # selected_patient_id 설정
    SLselected_gene = f'SLselected_gene <- "{SLselected_gene}"\n'
    r_script.insert(0, SLselected_gene)  # 스크립트 시작 부분에 삽입

    with open("C:/Users/hasn0737/Desktop/r_study/SL/temp_SL.R", "w") as file:
        file.writelines(r_script)

    #file_path = f"C:/Users/hasn0737/Desktop/r_study/SL/csv_file/{SLselected_gene}.csv"

    subprocess_result = subprocess.run(
        ["C:/Program Files/R/R-4.4.1/bin/Rscript", "C:/Users/hasn0737/Desktop/r_study/SL/temp_SL.R"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        cwd = "C:/Users/hasn0737/Desktop/r_study/SL/"
    )

    # R 스크립트 실행 결과 확인
    if subprocess_result.returncode == 0:
        print("R script executed successfully.")
        print("Output:", subprocess_result.stdout)
    else:
        print("Error occurred while running R script.")
        print("Stderr:", subprocess_result.stderr)

    return SLselected_gene

def SLselected_gene_input(request):

    context = {}
    if request.method == 'POST':

        SLselected_gene = request.POST.get('SLselected_gene')
        run_r_SL(SLselected_gene)
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


def r_script_process(request):
    try:
        # Celery 작업을 비동기적으로 실행
        task_result = run_r_script_task.apply_async()

        # 즉시 작업 ID 반환
        return JsonResponse({
            'status': 'pending',
            'task_id': task_result.id,
            'message': 'R script is being executed...'
        })

    except Exception as e:
        print(f"Error occurred: {str(e)}")
        return JsonResponse({
            'status': 'error',
            'message': 'An unexpected error occurred.',
            'error': str(e)
        })


# 작업 상태 조회 API
def r_script_process(request):
    try:
        # Celery 작업 실행
        task_result = run_r_script_task.apply_async()
        task_result.wait()  # 작업 완료까지 대기

        # 작업 상태 확인
        if task_result.state != 'SUCCESS':
            return JsonResponse({'status': 'error', 'message': 'Celery task did not complete successfully.', 'state': task_result.state})

        result = task_result.result
        if result and 'returncode' in result:
            if result['returncode'] == 0:
                return JsonResponse({'status': 'success', 'message': 'R script executed successfully.'})
            else:
                return JsonResponse({
                    'status': 'error',
                    'message': 'R script failed.',
                    'stdout': result['stdout'],
                    'stderr': result['stderr']
                })
        else:
            return JsonResponse({'status': 'error', 'message': 'No result returned from the task.'})

    except Exception as e:
        print(f"Error occurred: {str(e)}")
        return JsonResponse({
            'status': 'error',
            'message': 'An unexpected error occurred.',
            'error': str(e)
        })