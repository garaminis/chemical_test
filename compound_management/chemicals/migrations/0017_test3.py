# Generated by Django 4.1 on 2024-08-21 02:40

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("chemicals", "0016_test2"),
    ]

    operations = [
        migrations.CreateModel(
            name="Test3",
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
            options={
                "db_table": "test3",
            },
        ),
    ]
