# Generated by Django 5.0.1 on 2024-06-25 11:21

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('chemicals', '0016_rename_field_1_schrödingermodel_horalabs_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='chemical',
            name='MW',
            field=models.FloatField(blank=True, null=True),
        ),
    ]
