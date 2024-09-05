from django import forms
from django.contrib.auth.hashers import make_password
from django.contrib.auth.models import BaseUserManager, AbstractBaseUser, PermissionsMixin
from django.db import models
from django.utils import timezone
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

from users.models import User
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes.fields import GenericForeignKey, GenericRelation


class Chemical(models.Model):
    target = models.CharField(max_length=50)
    chem_id = models.CharField(max_length=50, primary_key=True)
    smiles = models.TextField()
    image = models.ImageField(upload_to='chemicals/', blank=True, null=True)
 #   image = models.ImageField(upload_to='chemical_images/', blank=True, null=True)
    MW = models.FloatField(null=True, blank=True)
    cLogP = models.FloatField(null=True, blank=True)
    TPSA = models.FloatField(blank=True, null=True)
    H_donors = models.IntegerField(blank=True, null=True)
    H_acceptors = models.IntegerField(blank=True, null=True)
    lipinski = models.BooleanField(default=False)
    user = models.CharField(max_length=200, null=True)
    datetime = models.DateTimeField(default=timezone.now, blank=False)

    def __str__(self):
        return f"{self.chem_id} - {self.smiles} - {self.target}"

    def save(self, *args, **kwargs):
        mol = Chem.MolFromSmiles(self.smiles)
        if mol:
            self.MW = round(Descriptors.MolWt(mol), 1)
            self.cLogP = round(Crippen.MolLogP(mol), 1)
            self.TPSA = round(Descriptors.TPSA(mol), 1)
            self.H_donors = Descriptors.NumHDonors(mol)
            self.H_acceptors = Descriptors.NumHAcceptors(mol)
            # 리핀스키의 규칙 체크
            self.lipinski = (self.MW < 500 and
                             self.cLogP < 5 and
                             self.H_donors <= 5 and
                             self.H_acceptors <= 10)
            # 분자 이미지 생성
            image_data = generate_image(self.smiles)
            # if image_data:
            #     self.image.save(f'{self.chem_id}.png', ContentFile(image_data), save=False)
        super().save(*args, **kwargs)

# 이미지 생성 함수
from rdkit.Chem import Draw
from io import BytesIO
from django.core.files.base import ContentFile
def generate_image(smiles, size=(1800, 1800)):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        return buffer.getvalue()
    return None

class Result(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    description = models.TextField()
    value = models.FloatField()
    timestamp = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.chemical} - {self.description}"

class Document(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE, related_name='documents')
    file = models.FileField(upload_to='uploads/')
    uploaded_at = models.DateTimeField(default=timezone.now)

class Pharmacokinetic(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    max_concentration = models.FloatField(null=True, blank=True)
    tmax = models.FloatField(null=True, blank=True)
    AUC = models.FloatField(null=True, blank=True)
    t_half = models.FloatField(null=True, blank=True)  # t1/2
    Vss = models.FloatField(null=True, blank=True)
    F = models.FloatField(null=True, blank=True)
    # BA = models.FloatField()  # Bioavailability
    CL = models.FloatField(null=True, blank=True)
    Route = models.CharField(max_length=200, null=True)
    user = models.CharField(max_length=200, null=True)
    files = GenericRelation('result_document')

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
    user = models.CharField(max_length=200, null=True)
    files = GenericRelation('result_document')

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date}'

class SchrödingerModel(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
#    QPlogPow = models.FloatField
    QPlogS = models.FloatField()
    QPlogHERG = models.FloatField()
    QPPCaco = models.FloatField()
    QPlogBB = models.FloatField()
    QPPMDCK = models.FloatField()
    Metab = models.FloatField()
    QPlogKhsa = models.FloatField()
    HOralAbs = models.FloatField()
    PerHOralAbs = models.FloatField()
    field_11 = models.FloatField()

    def __str__(self):
        return f'{self.chemical.chem_id} Schrödinger Model'

class LiverMicrosomalStability(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    mouse = models.FloatField()
    rat = models.FloatField()
    human = models.FloatField()
    user = models.CharField(max_length=200, null=True)
    files = GenericRelation('result_document')

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} Liver Microsomal Stability'

class CYPInhibition(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    cyp_1a2 = models.FloatField()
    cyp_2c9 = models.FloatField()
    cyp_2c19 = models.FloatField()
    cyp_2d6 = models.FloatField()
    cyp_3a4 = models.FloatField()
    user = models.CharField(max_length=200, null=True)
    files = GenericRelation('result_document')

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} CYP Inhibition'

class CCK_assay(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    user = models.CharField(max_length=200, null=True)
    cell = models.TextField(null=True, help_text='default')
    IC50 = models.FloatField(null=True, blank=True, help_text='check')
    Out = models.CharField(max_length=200,null=True, help_text='check')
    comment = models.TextField(null=True, help_text='default')
    files = GenericRelation('result_document')

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} CCK_assay'

class Western_blot(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    user = models.CharField(max_length=200, null=True)
    comment = models.TextField(null=True, help_text='default')
    files = GenericRelation('result_document')

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} Western_blot'

class Target_Inhibition(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    user = models.CharField(max_length=200, null=True)
    vitro_at_10 = models.FloatField(null=True, blank=True, help_text='check')
    Taret_IC50 = models.FloatField(null=True, blank=True, help_text='check')
    comment = models.TextField(null=True, help_text='default')
    files = GenericRelation('result_document')

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} Target_Inhibition'

class other_asssay(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    user = models.CharField(max_length=200, null=True)
    title = models.CharField(max_length=200, default="NA")
    comment = models.TextField(null=True)
    files = GenericRelation('result_document')

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} other_asssay'

class invtro_Image(models.Model):
    CCK_assay = models.ForeignKey(CCK_assay, related_name='images', on_delete=models.CASCADE,null=True)
    image = models.ImageField(upload_to='in_vitro/')
    Western_blot = models.ForeignKey(Western_blot, related_name='images', on_delete=models.CASCADE,null=True)
    Target_Inhibition = models.ForeignKey(Target_Inhibition,related_name='images',on_delete=models.CASCADE,null=True)
    other_asssay = models.ForeignKey(other_asssay,related_name='images',on_delete=models.CASCADE,null=True)

class in_vivo(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    start_date = models.DateField()
    end_date = models.DateField()
    user = models.CharField(max_length=200, null=True)
    cell = models.TextField(null=True)
    does = models.FloatField(null=True, blank=True)
    solvent =  models.TextField(null=True)
    inject_date= models.IntegerField(null=True)
    group = models.TextField(null=True)
    comment = models.TextField(null=True)
    category = models.CharField(max_length=200, null=True)
    files = GenericRelation('result_document')

    def __str__(self):
        return f'{self.chemical.chem_id} - in_vivo'

class Favorite(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    item = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        unique_together = ('user', 'item')



class FDA_result(models.Model):

    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    user = models.CharField(max_length=200, null=True)
    tmax = models.CharField(max_length=200, null=True)
    max_concentration = models.CharField(max_length=200, null=True)
    AUC = models.CharField(max_length=200, null=True)
    t_half = models.CharField(max_length=200, null=True)
    period = models.CharField(max_length=3, default='1')

    def __str__(self):
        return self.chemical

class result_document(models.Model):
    content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE)
    object_id = models.PositiveIntegerField()
    content_object = GenericForeignKey('content_type', 'object_id')

    file = models.FileField(upload_to='uploads/')
    uploaded_at = models.DateTimeField(default=timezone.now)

    def __str__(self):
        return self.file.name