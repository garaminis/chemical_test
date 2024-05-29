# Generated by Django 5.0.1 on 2024-05-28 23:41

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('chemicals', '0011_cytotoxicity'),
    ]

    operations = [
        migrations.CreateModel(
            name='SchrödingerModel',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('field_1', models.FloatField()),
                ('field_2', models.FloatField()),
                ('field_3', models.FloatField()),
                ('field_4', models.FloatField()),
                ('field_5', models.FloatField()),
                ('field_6', models.FloatField()),
                ('field_7', models.FloatField()),
                ('field_8', models.FloatField()),
                ('field_9', models.FloatField()),
                ('field_10', models.FloatField()),
                ('field_11', models.FloatField()),
                ('field_12', models.FloatField()),
                ('field_13', models.FloatField()),
                ('field_14', models.FloatField()),
                ('field_15', models.FloatField()),
                ('field_16', models.FloatField()),
                ('field_17', models.FloatField()),
                ('field_18', models.FloatField()),
                ('field_19', models.FloatField()),
                ('field_20', models.FloatField()),
                ('chemical', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemicals.chemical')),
            ],
        ),
    ]
