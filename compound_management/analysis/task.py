from celery import shared_task
import subprocess


@shared_task
def run_r_script_task():
   result = subprocess.run(
        ["C:/Program Files/R/R-4.4.1/bin/Rscript", "C:/Users/hasn0737/Desktop/r_study/SL/SL_start.R"],
        capture_output=True,
        text=True,
        encoding="utf-8"
    )
   return {
       'stdout': result.stdout,
       'stderr': result.stderr,
       'returncode': result.returncode
   }

@shared_task
def test_task():
    return 'Task completed!'