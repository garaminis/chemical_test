# Generated by Django 5.0.1 on 2024-05-26 12:29

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('chemicals', '0002_alter_chemical_chem_id_alter_chemical_str_img'),
    ]

    operations = [
        migrations.AlterField(
            model_name='chemical',
            name='cLogP',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='chemical',
            name='smiles',
            field=models.TextField(),
        ),
        migrations.AlterField(
            model_name='chemical',
            name='str_img',
            field=models.ImageField(blank=True, null=True, upload_to='images/'),
        ),
    ]
