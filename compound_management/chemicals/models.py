<<<<<<< HEAD
from django import forms
=======
>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf
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
<<<<<<< HEAD
=======

>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf

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
<<<<<<< HEAD
    max_concentration = models.FloatField(null=True, blank=True)
=======
    cmax = models.FloatField(null=True, blank=True)
>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf
    tmax = models.FloatField(null=True, blank=True)
    AUC = models.FloatField(null=True, blank=True)
    t_half = models.FloatField(null=True, blank=True)  # t1/2
    Vss = models.FloatField(null=True, blank=True)
    F = models.FloatField(null=True, blank=True)
    # BA = models.FloatField()  # Bioavailability
    CL = models.FloatField(null=True, blank=True)
    Route = models.CharField(max_length=200, null=True)
    user = models.CharField(max_length=200, null=True)
<<<<<<< HEAD
    files = GenericRelation('result_document')
=======
>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf

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
<<<<<<< HEAD
    files = GenericRelation('result_document')
=======
>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf

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
<<<<<<< HEAD
    files = GenericRelation('result_document')
=======
>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf

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
<<<<<<< HEAD
    files = GenericRelation('result_document')
=======
>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} CYP Inhibition'

class CCK_assay(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    user = models.CharField(max_length=200, null=True)
<<<<<<< HEAD
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
=======
    cell = models.TextField(null=True)
    IC50 = models.FloatField(null=True, blank=True)
    Out = models.CharField(max_length=200,null=True)
    comment = models.TextField(null=True)

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} CCK_assay'

class Western_blot(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    user = models.CharField(max_length=200, null=True)
    comment = models.TextField(null=True)

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} Western_blot'

class Target_Inhibition(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    user = models.CharField(max_length=200, null=True)
    vitro_at_10 = models.FloatField(null=True, blank=True)
    Taret_IC50 = models.FloatField(null=True, blank=True)
    comment = models.TextField(null=True)

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} Target_Inhibition'

class other_asssay(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    user = models.CharField(max_length=200, null=True)
    title = models.CharField(max_length=200, default="NA")
    comment = models.TextField(null=True)

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} other_asssay'


class invtro_Image(models.Model):
    cck_assay = models.ForeignKey(CCK_assay, related_name='images', on_delete=models.CASCADE,null=True)
    image = models.ImageField(upload_to='in_vitro/')
    wb = models.ForeignKey(Western_blot, related_name='images', on_delete=models.CASCADE,null=True)
    in_target = models.ForeignKey(Target_Inhibition,related_name='images',on_delete=models.CASCADE,null=True)
    others = models.ForeignKey(other_asssay,related_name='images',on_delete=models.CASCADE,null=True)

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

    def __str__(self):
        return f'{self.chemical.chem_id} - in_vivo'
class UserManager(BaseUserManager):
    def create_user(self, email, userID, password=None, **extra_fields):

        if not email:
            raise ValueError('The Email field must be set')
        if not userID:
            raise ValueError('The userID field must be set')
        email = self.normalize_email(email)
        user = self.model(email=email, userID=userID, **extra_fields)
        user.set_password(password)
        user.save(using=self._db)
        return user

    def create_superuser(self, email, userID, password=None, **extra_fields):

        extra_fields.setdefault('is_staff', True)
        extra_fields.setdefault('is_superuser', True)

        if extra_fields.get('is_staff') is not True:
            raise ValueError('Superuser must have is_staff=True.')
        if extra_fields.get('is_superuser') is not True:
            raise ValueError('Superuser must have is_superuser=True.')
        return self.create_user(email, userID, password, **extra_fields)

class User(AbstractBaseUser,PermissionsMixin):
    id = models.AutoField(primary_key=True)
    userID = models.CharField(default='', max_length=100, null=False, blank=False, unique=True)
    email = models.EmailField(default='', max_length=100, null=False, blank=False, unique=True)
    name = models.CharField(default='', max_length=100, null=False, blank=False)
    roll = models.CharField(default='', max_length=100, blank=False)
    group = models.CharField(default='', max_length=100, blank=False)
    date_joined = models.DateTimeField(default=timezone.now, blank=False)
    is_staff = models.BooleanField(default=False)
    is_active = models.BooleanField(default=True)
    is_superuser = models.BooleanField(default=False)

    objects = UserManager()
    USERNAME_FIELD = 'userID'
    REQUIRED_FIELDS = ['email', 'name']

    def __str__(self): # 객체의 문자열 표현을 반환하는 메서드
        return self.userID

    def has_perm(self, perm, obj=None): # 슈퍼유저는 모든 권한을 가짐.
        return self.is_superuser

    # def has_module_perms(self, app_label): # 모든 app에서 슈퍼유저의 권한 부여
    #     return self.is_superuser
>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf
