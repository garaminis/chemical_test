# Generated by Django 4.1 on 2024-08-30 05:33

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):
    dependencies = [
        ("chemicals", "0003_fda_result_file_alter_fda_result_period"),
    ]

    operations = [
        migrations.RemoveField(
            model_name="fda_result",
            name="file",
        ),
        migrations.CreateModel(
            name="Document",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("file", models.FileField(upload_to="uploads/")),
                ("uploaded_at", models.DateTimeField(auto_now_add=True)),
                (
                    "chemical",
                    models.ForeignKey(
                        on_delete=django.db.models.deletion.CASCADE,
                        to="chemicals.chemical",
                    ),
                ),
            ],
        ),
    ]
