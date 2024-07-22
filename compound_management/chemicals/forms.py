from django import forms
from .models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel,LiverMicrosomalStability, CYPInhibition ,User


class DateInput(forms.DateInput):
    input_type = 'date'


class ChemicalForm(forms.ModelForm):
    class Meta:
        model = Chemical
        fields = ['chem_id', 'smiles', 'image', 'MW', 'cLogP', 'TPSA','H_donors', 'H_acceptors', 'lipinski' ]

class ChemicalUploadForm(forms.Form):
    file = forms.FileField()

class PharmacokineticForm(forms.ModelForm):
    class Meta:
        model = Pharmacokinetic
        fields = ['date', 'cmax', 'tmax', 'AUC', 't_half', 'Vss', 'Vd', 'BA']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class CytotoxicityForm(forms.ModelForm):
    class Meta:
        model = Cytotoxicity
        fields = ['date', 'VERO', 'HFL_1', 'L929', 'NIH_3T3', 'CHO_K1']
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
        fields = ['date', 'mouse', 'rat', 'human']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class CYPInhibitionForm(forms.ModelForm):
    class Meta:
        model = CYPInhibition
        fields = ['date', 'cyp_1a2', 'cyp_2c9', 'cyp_2c19', 'cyp_2d6', 'cyp_3a4']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class UserForm(forms.ModelForm):
    # 패스워드 필드를 정의하며, 입력 시 비밀번호 입력 필드를 사용합니다.
    password = forms.CharField(widget=forms.PasswordInput())
    password2 = forms.CharField(widget=forms.PasswordInput())
    email = forms.EmailField()

    class Meta:
        # User 모델을 기반으로 폼을 정의합니다.
        model = User
        # 폼에서 사용할 모델 필드를 지정합니다.
        fields = ['userID', 'email', 'name', 'password','password2', 'roll', 'group']

    def clean_userID(self):
        userID = self.cleaned_data.get('userID')
        if User.objects.filter(userID=userID).exists():
            raise forms.ValidationError("중복된 아이디 입니다.")
        return userID

    def clean_password(self):
        password = self.cleaned_data.get('password')
        if len(password) < 8:
            raise forms.ValidationError("Password must be at least 8 characters long.")
        return password

    def clean_password2(self):
        # 두 비밀번호 입력 일치 확인
        password = self.cleaned_data.get("password")
        password2 = self.cleaned_data.get("password2")
        if password and password2 and password != password2:
            raise forms.ValidationError("Passwords don't match")
        return password
    #
    # def save(self, commit=True):
    #     # 기본 save 메서드를 호출하여 User 인스턴스를 가져옵니다.
    #     user = super().save(commit=False)
    #     # 입력받은 패스워드를 해시화하여 User 인스턴스에 설정합니다.
    #     user.set_password(self.cleaned_data['password'])
    #     # commit이 True일 경우, 데이터베이스에 저장합니다.
    #     if commit:
    #         user.save()
    #     # 저장된 User 인스턴스를 반환합니다.
    #     return user