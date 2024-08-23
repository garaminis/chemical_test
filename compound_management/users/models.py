from django.contrib.auth.base_user import BaseUserManager, AbstractBaseUser
from django.contrib.auth.models import PermissionsMixin
from django.db import models
from django.utils import timezone

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


class DatabaseList(models.Model):
    name = models.CharField(max_length=100, unique=True)
    group = models.CharField(max_length=100, null=True)
    title = models.CharField(max_length=100 ,unique=True)

    def __str__(self):
        return self.name

