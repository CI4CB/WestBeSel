"""
@author: jibril mammeri | heig-vd
"""

from django.core.management.base import BaseCommand, CommandError
from prediction.models import *
from datetime import datetime
import os

class Command(BaseCommand):

    def handle(self, *args, **options):
        proc=Process.objects.filter(exp_date__lte=datetime.now())
        for i in proc:
            name=i.name+'_'+str(i.id)
            print name
            path1='/code/westbesel/static/prediction/output/'+name
            path2='/code/westbesel/static/prediction/RESULTS_test/'+name+"_aaCS*"
            com='rm -rf '
            os.system(com+path1)
            os.system(com+path2)
            i.delete()
