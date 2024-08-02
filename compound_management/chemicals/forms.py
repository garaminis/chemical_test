import re

from django import forms
from django.forms import modelformset_factory

from .models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel, LiverMicrosomalStability, CYPInhibition, \
    User, CCK_assay, CCK_Image


# from compound_management.chemicals.models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel, LiverMicrosomalStability, CYPInhibition, User


class DateInput(forms.DateInput):
    input_type = 'date'

class ChemicalForm(forms.ModelForm):
    class Meta:
        model = Chemical
        fields = ['chem_id', 'smiles', 'image', 'MW', 'cLogP', 'TPSA','H_donors', 'H_acceptors', 'lipinski','user']

class ChemicalUploadForm(forms.Form):
    file = forms.FileField()

class PharmacokineticForm(forms.ModelForm):
    class Meta:
        model = Pharmacokinetic
        fields = ['date', 'cmax', 'tmax', 'AUC', 't_half', 'Vss', 'F', 'CL', 'Route','user']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class CytotoxicityForm(forms.ModelForm):
    class Meta:
        model = Cytotoxicity
        fields = ['date', 'VERO', 'HFL_1', 'L929', 'NIH_3T3', 'CHO_K1','user']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class SchrödingerModelForm(forms.ModelForm):
    class Meta:
        model = SchrödingerModel
        fields = ['QPlogS', 'QPlogHERG', 'QPPCaco', 'QPlogBB', 'QPPMDCK', 'Metab', 'QPlogKhsa', 'HOralAbs', 'PerHOralAbs', 'field_11']

class SchrödingerModelUploadForm(forms.Form):
    file = forms.FileField()

class LiverMicrosomalStabilityForm(forms.ModelForm):
    class Meta:
        model = LiverMicrosomalStability
        fields = ['date', 'mouse', 'rat', 'human','user']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class CYPInhibitionForm(forms.ModelForm):
    class Meta:
        model = CYPInhibition
        fields = ['date', 'cyp_1a2', 'cyp_2c9', 'cyp_2c19', 'cyp_2d6', 'cyp_3a4','user']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class cckForm(forms.ModelForm):
    class Meta:
        model = CCK_assay
        fields = ['date', 'user',  'cell', 'IC50' , 'Out', 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }
class cckimgForm(forms.ModelForm):
    class Meta:
        model = CCK_Image
        fields = ['image']
        widgets = {
            'image': forms.ClearableFileInput(attrs={'multiple': True}),
        }

class UserForm(forms.ModelForm):
    # 패스워드 필드를 정의하며, 입력 시 비밀번호 입력 필드를 사용합니다.
    password = forms.CharField(widget=forms.PasswordInput)
    password2 = forms.CharField(widget=forms.PasswordInput())

    class Meta:
        model = User
        fields = ['userID', 'email', 'name', 'password', 'password2', 'roll', 'group']

    def clean_password(self):
        password = self.cleaned_data.get('password')
        pattern = r'^(?=.*[A-Za-z])(?=.*\d)(?=.*[@$!%*#?&])[A-Za-z\d@$!%*#?&]{8}$'
        if not re.match(pattern, password ):
            raise forms.ValidationError("문자, 숫자, 특수문자를 포함한 8자리 이상이어야 합니다.")
        return password

    def clean_userID(self):
        userID = self.cleaned_data.get('userID')
        if User.objects.filter(userID=userID).exists():
            raise forms.ValidationError("중복된 아이디입니다.")
        if len(userID) < 4:
            print(userID)
            raise forms.ValidationError("아이디는 4자 이상입니다.")
        return userID

    def clean_email(self):
        email = self.cleaned_data.get('email')
        # 이메일 형식 검증
        if not self.is_valid_email_format(email):
            raise forms.ValidationError('유효하지 않은 이메일 주소 형식입니다.')
        return email
    def is_valid_email_format(self, email):
        # 이메일 형식 검증을 위한 정규 표현식
        email_regex = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
        return re.match(email_regex, email) is not None

    def save(self, commit=True):
        user = super().save(commit=False) #save메서드 호출,커밋 하지 않은 상태로 user생성
        user.set_password(self.cleaned_data['password']) # 비밀번호 가져와서 set_password메서드 사용하여 해시화
        if commit: #커밋이 true이면 저장
            user.save()
        return user


