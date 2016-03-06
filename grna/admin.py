from django.contrib import admin

# Register your models here.
from .models import *

admin.site.register(Species)
admin.site.register(Target)
admin.site.register(Nuclease)
admin.site.register(PAM)
admin.site.register(GuideRNA)
