from django import forms
<<<<<<< HEAD
from django.forms import ClearableFileInput, formset_factory
from .models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel, LiverMicrosomalStability, CYPInhibition, \
    CCK_assay, invtro_Image, Western_blot, Target_Inhibition, other_asssay, in_vivo, FDA_result, Document
=======
from django.forms import modelformset_factory, ClearableFileInput

from .models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel, LiverMicrosomalStability, CYPInhibition, \
    User, CCK_assay, invtro_Image, Western_blot, Target_Inhibition, other_asssay, in_vivo


# from compound_management.chemicals.models import Chemical, Pharmacokinetic, Cytotoxicity, SchrödingerModel, LiverMicrosomalStability, CYPInhibition, User
>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf


class DateInput(forms.DateInput):
    input_type = 'date'

class ChemicalForm(forms.ModelForm):
    class Meta:
        model = Chemical
        fields = ['chem_id', 'smiles', 'image', 'MW', 'cLogP', 'TPSA' , 'H_donors', 'H_acceptors', 'lipinski','user']
<<<<<<< HEAD

class DocumentForm(forms.ModelForm):
    class Meta:
        model = Document
        fields = ['file']
=======
>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf

class ChemicalUploadForm(forms.Form):
    file = forms.FileField()

class PharmacokineticForm(forms.ModelForm):
    class Meta:
        model = Pharmacokinetic
<<<<<<< HEAD
        fields = ['date', 'max_concentration', 'tmax', 'AUC', 't_half', 'Vss', 'F', 'CL', 'Route','user']
=======
        fields = ['date', 'cmax', 'tmax', 'AUC', 't_half', 'Vss', 'F', 'CL', 'Route','user']
>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf
        widgets = {
            'date': DateInput(attrs={'type': 'date'}),  # DateInput 위젯 사용
        }

class CytotoxicityForm(forms.ModelForm):
    class Meta:
        model = Cytotoxicity
        fields = ['date', 'VERO', 'HFL_1', 'L929', 'NIH_3T3', 'CHO_K1', 'user']
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
        fields = ['date', 'mouse', 'rat', 'human' , 'user']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class CYPInhibitionForm(forms.ModelForm):
    class Meta:
        model = CYPInhibition
        fields = ['date', 'cyp_1a2', 'cyp_2c9', 'cyp_2c19', 'cyp_2d6', 'cyp_3a4' , 'user']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class cckForm(forms.ModelForm):
<<<<<<< HEAD
    class Meta:
        model = CCK_assay
        fields = ['date', 'user',  'cell', 'IC50' , 'Out' , 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }
class wbForm(forms.ModelForm):
    class Meta:
        model = Western_blot
        fields = ['date', 'user' , 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }
=======
    class Meta:
        model = CCK_assay
        fields = ['date', 'user',  'cell', 'IC50' , 'Out' , 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }
class wbForm(forms.ModelForm):
    class Meta:
        model = Western_blot
        fields = ['date', 'user' , 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class intargetForm(forms.ModelForm):
    class Meta:
        model = Target_Inhibition
        fields = ['date', 'user', 'vitro_at_10' , 'Taret_IC50' , 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class otherForm(forms.ModelForm):
    class Meta:
        model = other_asssay
        fields = ['date', 'user', 'title' , 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

# class invtroimgForm(forms.ModelForm):
#     images = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True}))
#     class Meta:
#         model = invtro_Image
#         fields = ['image']

class in_vivoForm(forms.ModelForm):
    class Meta:
        model = in_vivo
        fields = ['user', 'start_date', 'end_date', 'cell', 'does', 'solvent', 'inject_date', 'group', 'comment', 'category']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }
class CustomClearableFileInput(ClearableFileInput):
    allow_multiple_selected = True

    def __init__(self, attrs=None):
        super().__init__(attrs)
        if attrs is None:
            attrs = {}
        attrs.update({'multiple': True})

    def value_from_datadict(self, data, files, name):
        if self.allow_multiple_selected:
            return files.getlist(name)
        return files.get(name)
class InvtroimgForm(forms.ModelForm):
    images = forms.ImageField(widget=CustomClearableFileInput)
    class Meta:
        model = invtro_Image # 여기에 해당 모델 이름을 입력하세요
        fields = ['images']

class UserForm(forms.ModelForm):
    # 패스워드 필드를 정의하며, 입력 시 비밀번호 입력 필드를 사용합니다.
    password = forms.CharField(widget=forms.PasswordInput)
    password2 = forms.CharField(widget=forms.PasswordInput())

    class Meta:
        model = User
        fields = ['userID', 'email', 'name', 'password', 'password2', 'roll', 'group']
>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf

class intargetForm(forms.ModelForm):
    class Meta:
        model = Target_Inhibition
        fields = ['date', 'user', 'vitro_at_10' , 'Taret_IC50' , 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

class otherForm(forms.ModelForm):
    class Meta:
        model = other_asssay
        fields = ['date', 'user', 'title' , 'comment']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }

# class invtroimgForm(forms.ModelForm):
#     images = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True}))
#     class Meta:
#         model = invtro_Image
#         fields = ['image']

<<<<<<< HEAD
class in_vivoForm(forms.ModelForm):
    class Meta:
        model = in_vivo
        fields = ['user', 'start_date', 'end_date', 'cell', 'does', 'solvent', 'inject_date', 'group', 'comment', 'category']
        widgets = {
            'date': DateInput(),  # DateInput 위젯 사용
        }
class CustomClearableFileInput(ClearableFileInput):
    allow_multiple_selected = True
=======
    def save(self, commit=True):
        user = super().save(commit=False) #save메서드 호출,커밋 하지 않은 상태로 user생성
        user.set_password(self.cleaned_data['password']) # 비밀번호 가져와서 set_password메서드 사용하여 해시화
        if commit: #커밋이 true이면 저장
            user.save()
        return user

>>>>>>> 323fb36a92ed679bfac130e81025e71786d36baf

    def __init__(self, attrs=None):
        super().__init__(attrs)
        if attrs is None:
            attrs = {}
        attrs.update({'multiple': True})

    def value_from_datadict(self, data, files, name):
        if self.allow_multiple_selected:
            return files.getlist(name)
        return files.get(name)

class InvtroimgForm(forms.ModelForm):
    images = forms.ImageField(widget=CustomClearableFileInput)
    class Meta:
        model = invtro_Image
        fields = ['images']

class FDA_Form(forms.ModelForm):
    class Meta:
        model = FDA_result
        fields = ['user','tmax','max_concentration','AUC','t_half','period']

class FDA_UploadForm(forms.Form):
    file = forms.FileField()
