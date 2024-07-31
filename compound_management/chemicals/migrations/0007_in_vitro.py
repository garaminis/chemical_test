# Generated by Django 4.1 on 2024-07-31 04:02

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):
    dependencies = [
        ("chemicals", "0006_chemical_user_pharmacokinetic_user"),
    ]

    operations = [
        migrations.CreateModel(
            name="In_vitro",
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
                ("date", models.DateField()),
                ("user", models.CharField(max_length=200, null=True)),
                ("assay", models.CharField(max_length=200, null=True)),
                ("cell", models.TextField()),
                ("IC50", models.FloatField()),
                ("vitro_at_10", models.FloatField()),
                ("Taret_IC50", models.FloatField()),
                ("Out", models.BooleanField()),
                (
                    "image",
                    models.ImageField(blank=True, null=True, upload_to="in_vitro/"),
                ),
                ("comment", models.TextField()),
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
