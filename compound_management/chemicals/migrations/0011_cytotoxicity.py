# Generated by Django 5.0.1 on 2024-05-28 14:50

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('chemicals', '0010_pharmacokinetic'),
    ]

    operations = [
        migrations.CreateModel(
            name='Cytotoxicity',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('date', models.DateField()),
                ('VERO', models.FloatField()),
                ('HFL_1', models.FloatField()),
                ('L929', models.FloatField()),
                ('NIH_3T3', models.FloatField()),
                ('CHO_K1', models.FloatField()),
                ('chemical', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chemicals.chemical')),
            ],
        ),
    ]
