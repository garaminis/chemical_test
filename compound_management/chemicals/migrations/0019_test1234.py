# Generated by Django 4.1 on 2024-08-21 07:43

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("chemicals", "0018_delete_other123_delete_test2_delete_test3_and_more"),
    ]

    operations = [
        migrations.CreateModel(
            name="Test1234",
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
                "db_table": "test1234",
            },
        ),
    ]
