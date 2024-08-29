import os

import django
from django.contrib import messages
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.shortcuts import render, redirect, get_object_or_404
from django.views.decorators.http import require_POST

from .forms import UserForm



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


# @login_required
# @require_POST
# def toggle_favorite(request, chem_id):
#     print(f"Toggle favorite called with chem_id: {chem_id}")  # 이 줄을 추가
#     chemical = get_object_or_404(Chemical, chem_id=chem_id)
#     favorite, created = Favorite.objects.get_or_create(user=request.user, item=chemical)
#
#     if not created:
#         favorite.delete()
#         is_favorite = False
#     else:
#         is_favorite = True
#
#     return JsonResponse({'success': True, 'is_favorite': is_favorite})


