from django.contrib import admin
from django.urls import path, include
from django.conf import settings
from django.conf.urls.static import static


urlpatterns = [
    path('admin/', admin.site.urls),
    path('', include('chemicals.urls')),
    path('users/', include('users.urls')),
    path('data/', include('data.urls')), # data 앱의 URL 패턴을 인식하고 처리
]
if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
