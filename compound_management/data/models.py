from django.contrib.contenttypes.fields import GenericRelation
from django.db import models
from chemicals.models import Chemical, result_document



# Create your models here.

class other_assay1 (models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    user = models.CharField(max_length=200, null=True)
    title = models.CharField(max_length=200, default="NA")
    comment = models.TextField(null=True,help_text='default')
    files = GenericRelation(result_document)

    class Meta:
        app_label = 'data'


class other_assay2 (models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    user = models.CharField(max_length=200, null=True)
    title = models.CharField(max_length=200, default="NA",help_text='default')
    HFL_1 = models.FloatField(help_text='check')
    L929 = models.FloatField(help_text='check')
    files = GenericRelation(result_document)
    class Meta:
        app_label = 'data'