# Generated by Django 4.1 on 2024-08-13 07:08

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):
    dependencies = [
        ("chemicals", "0025_chemical_datetime"),
    ]

    operations = [
        migrations.CreateModel(
            name="in_vivo",
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
                ("start_date", models.DateField()),
                ("end_date", models.DateField()),
                ("user", models.CharField(max_length=200, null=True)),
                ("cell", models.TextField(null=True)),
                ("does", models.FloatField(blank=True, null=True)),
                ("solvent", models.TextField(null=True)),
                ("inject_date", models.IntegerField(null=True)),
                ("group", models.TextField(null=True)),
                ("comment", models.TextField(null=True)),
                ("category", models.CharField(max_length=200, null=True)),
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
