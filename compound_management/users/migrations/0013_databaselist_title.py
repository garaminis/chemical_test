# Generated by Django 4.1 on 2024-08-23 00:21

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("users", "0012_create_teat2"),
    ]

    operations = [
        migrations.AddField(
            model_name="databaselist",
            name="title",
            field=models.CharField(max_length=100, null=True),
        ),
    ]
