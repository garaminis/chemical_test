# Generated by Django 4.1 on 2024-08-21 08:15

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("chemicals", "0021_delete_test1234_delete_test22"),
    ]

    operations = [
        migrations.CreateModel(
            name="Test123456",
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
                ("name", models.CharField(max_length=255)),
            ],
        ),
    ]
