import os
from celery import Celery

# compound_management.settings을 Celery에서 사용하도록 설정
os.environ.setdefault(
    'DJANGO_SETTINGS_MODULE',
    'compound_management.settings'
)
os.environ.setdefault('FORKED_BY_MULTIPROCESSING', '1') # 호환성 충돌 방지!

# Celery 생성
app = Celery('compound_management')

# Celery 관련 설정을 읽어옴
app.config_from_object(
    'django.conf:settings',
    namespace='CELERY')

# 자동 Task 등록함
app.autodiscover_tasks()
