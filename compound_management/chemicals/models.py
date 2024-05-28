from django.db import models

# Create your models here.
from django.contrib.auth.models import User

# User 모델을 확장하거나 기본 User 모델을 사용할 수 있습니다.
class Chemical(models.Model):
    target = models.CharField(max_length=50)
    chem_id = models.CharField(max_length=50, primary_key=True)
    smiles = models.TextField()
    image = models.ImageField(upload_to='chemicals/', blank=True, null=True)
 #   image = models.ImageField(upload_to='chemical_images/', blank=True, null=True)
    MW = models.FloatField()
    cLogP = models.FloatField(null=True, blank=True)


    def __str__(self):
        return f"{self.chem_id} - {self.smiles} - {self.target}"

class Result(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    description = models.TextField()
    value = models.FloatField()
    timestamp = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.chemical} - {self.description}"

class Pharmacokinetic(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    cmax = models.FloatField()
    tmax = models.FloatField()
    AUC = models.FloatField()
    t_half = models.FloatField()  # t1/2
    Vss = models.FloatField()
    Vd = models.FloatField()
    BA = models.FloatField()  # Bioavailability

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date}'

class Cytotoxicity(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    VERO = models.FloatField()
    HFL_1 = models.FloatField()
    L929 = models.FloatField()
    NIH_3T3 = models.FloatField()
    CHO_K1 = models.FloatField()

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date}'