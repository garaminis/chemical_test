import os
from celery import Celery

# compound_management.settings을 Celery에서 사용하도록 설정
os.environ.setdefault(
    'DJANGO_SETTINGS_MODULE',
    'compound_management.settings'
)
os.environ.setdefault(
    'FORKED_BY_MULTIPROCESSING',
    '1'
) # 호환성 충돌 방지!

# Celery 생성
app = Celery('compound_management')

# Celery 관련 설정을 읽어옴
app.config_from_object(
    'django.conf:settings',
    namespace='CELERY'
)
app.conf.update(
    worker_max_tasks_per_child=1,  # 각 워커가 하나의 태스크를 처리한 후 재시작되도록 설정
    worker_force_execv=True,       # execv를 강제 실행
)

# 자동 Task 등록함
app.autodiscover_tasks()

broker_connection_retry_on_startup = True
broker_connection_retry = True  # 현재 버전에서는 여전히 필요
broker_connection_retry_on_startup = True  # 향후 호환성을 위해 추가