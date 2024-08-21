import os

import django
from django.contrib import messages, admin
from django.contrib.admin.views.decorators import staff_member_required
from django.contrib.auth import authenticate, login, logout
from django.core.management import call_command
from django.db import models
from django.shortcuts import render, redirect
from django.utils.text import capfirst
from .forms import UserForm
from .models import DatabaseList
from django.db import connections
# from Chem.compound_management.users.forms import UserForm, TableForm, ColumnFormSet

def register_view(request):
    if request.method == 'POST':
        form = UserForm(request.POST)
        if form.is_valid():
            form.save()
            return redirect('login')
    else:
        form = UserForm()
    return render(request, 'chemicals/register.html', {'form': form}) # 오류 렌더링

def login_view(request):
    if request.method == 'POST':
        userID = request.POST['userID']     #userID = request.POST.get('userID')도 가능
        password = request.POST['password']
        user = authenticate(request, username=userID, password=password)
        print(user)
        if user is not None:
            login(request, user)
            return redirect('home')
        else:
            messages.error(request, '로그인에 실패하였습니다.')
    return render(request, 'chemicals/login.html')

def logout_view(request):
    logout(request)
    return redirect('login')

