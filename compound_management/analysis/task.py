import json
import time
from celery import shared_task
import subprocess
import pyRserve
import socket
from django.core.cache import cache
import logging

logger = logging.getLogger(__name__)

@shared_task()
def run_r_script_task():
    try:
        result = subprocess.run(
            ["C:/Program Files/R/R-4.4.1/bin/Rscript", "C:/Users/hasn0737/Desktop/r_study/SL/SL_start.R"],  # R 스크립트 경로
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8"
        )
        print(result.stdout)  # R 스크립트 출력 확인

    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e.stderr}")