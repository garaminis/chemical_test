from django.contrib.auth.hashers import make_password
from django.contrib.auth.models import BaseUserManager, AbstractBaseUser, PermissionsMixin
from django.db import models
from django.utils import timezone
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

# Create your models here.
# User 모델을 확장하거나 기본 User 모델을 사용할 수 있습니다.
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
            if image_data:
                self.image.save(f'{self.chem_id}.png', ContentFile(image_data), save=False)
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

class Pharmacokinetic(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    cmax = models.FloatField(null=True, blank=True)
    tmax = models.FloatField(null=True, blank=True)
    AUC = models.FloatField(null=True, blank=True)
    t_half = models.FloatField(null=True, blank=True)  # t1/2
    Vss = models.FloatField(null=True, blank=True)
    F = models.FloatField(null=True, blank=True)
    # BA = models.FloatField()  # Bioavailability
    CL = models.FloatField(null=True, blank=True)
    Route = models.CharField(max_length=200, null=True)
    user = models.CharField(max_length=200, null=True)

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

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} CYP Inhibition'

class In_vitro(models.Model):
    chemical = models.ForeignKey(Chemical, on_delete=models.CASCADE)
    date = models.DateField()
    user = models.CharField(max_length=200, null=True)
    assay = models.CharField(max_length=200, null=True)
    cell = models.TextField(null=True)
    IC50 = models.FloatField(null=True, blank=True)
    vitro_at_10 = models.FloatField(null=True, blank=True)
    Taret_IC50 = models.FloatField(null=True, blank=True)
    Out = models.CharField(max_length=200,null=True)
    image = models.ImageField(upload_to='in_vitro/', blank=True, null=True)
    comment = models.TextField(null=True)
    title = models.CharField(max_length=200,default="NA")

    def __str__(self):
        return f'{self.chemical.chem_id} - {self.date} In_vitro'

# class In_vitro_Image(models.Model):
#     in_vitro_id = models.ForeignKey(In_vitro, related_name='images', on_delete=models.CASCADE)
#     image = models.ImageField(upload_to='in_vitro/')
#
#     def __str__(self):
#         return f"Image for {self.chemical.name}"

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