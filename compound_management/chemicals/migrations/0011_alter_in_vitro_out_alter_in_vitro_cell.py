# Generated by Django 4.1 on 2024-07-31 06:00

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("chemicals", "0010_in_vitro_title"),
    ]

    operations = [
        migrations.AlterField(
            model_name="in_vitro",
            name="Out",
            field=models.BooleanField(default="NA", null=True),
        ),
        migrations.AlterField(
            model_name="in_vitro",
            name="cell",
            field=models.TextField(null=True),
        ),
    ]
