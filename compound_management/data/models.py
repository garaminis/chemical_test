from django.db import models

# Create your models here.


class Test0(models.Model):
    name = models.CharField(max_length=255)
    class Meta:
        db_table = 'test0'
        app_label = 'chemicals'
