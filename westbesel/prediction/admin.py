# -*- coding: utf-8 -*-
"""
@author: jibril mammeri | heig-vd
"""

from __future__ import unicode_literals

from django.contrib import admin
from .models import *
# Register your models here.


class ProcessAdmin(admin.ModelAdmin):
    list_display = ('id','name','step','_typ','user','mail','threshold','stretch','cystein','repeats','over',"fail",'date','exp_date')
    def _typ(self, obj):
        return "http://westbesel.iict.ch/results/"+obj.hash_id
admin.site.register(Process, ProcessAdmin)
