import json

from celery import shared_task
import subprocess
import pyRserve
from django.core.cache import cache
import logging
from rpy2.robjects import pandas2ri
from rpy2 import robjects

logger = logging.getLogger(__name__)

@shared_task
def run_r_script_task():
    try:
        # Rserve에 연결
        conn = pyRserve.connect(host='localhost', port=6311)

        # 첫 번째 R 스크립트 실행 (라이브러리 및 데이터 로드)
        conn.voidEval("source('C:/Users/hasn0737/Desktop/r_study/SL/SL_start.R')")

        # 필요한 데이터 변수 저장 (예: CSV 파일 처리 결과를 'result_data' 변수로 저장)
        result_data = conn.eval('result_data')

        # 캐시에 저장
        cache.set('r_script_result', result_data, timeout=3600)

        conn.close()  # 세션 닫기

        return {
            'stdout': result_data,
            'stderr': '',
            'returncode': 0
        }

    except Exception as e:
        logger.error(f"Error occurred while running R script: {str(e)}")
        return {
            'stdout': '',
            'stderr': str(e),
            'returncode': 1
        }