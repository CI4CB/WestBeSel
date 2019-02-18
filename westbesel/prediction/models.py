# -*- coding: utf-8 -*-
"""
@author: jibril mammeri | heig-vd
"""

from __future__ import unicode_literals

from django.db import models
from django.contrib.auth.models import User
from datetime import datetime

class Process(models.Model):
    name = models.CharField(max_length = 200)
    organism = models.CharField(max_length = 200,null = True, blank = True)
    protein = models.CharField(max_length = 200,null = True, blank = True)
    uniprot = models.CharField(max_length = 200,null = True, blank = True)
    hash_id = models.CharField(max_length = 200)
    step = models.CharField(max_length = 200)
    mail = models.CharField(max_length = 200,null = True, blank = True)
    min_ep_length = models.IntegerField(null = True, blank = True)
    max_ep_length = models.IntegerField(null = True, blank = True)
    min_ter_length = models.IntegerField(null = True, blank = True)
    max_ter_length = models.IntegerField(null = True, blank = True)
    threshold = models.FloatField(null = True, blank = True)
    stretch = models.IntegerField(null = True, blank = True)
    repeats = models.BooleanField(default = True, blank = True)
    cystein = models.BooleanField(default = True, blank = True)
    over = models.BooleanField(default = True, blank = True)
    date = models.DateTimeField(default = datetime.now, blank = True)
    exp_date = models.DateTimeField(null = True, blank = True)
    refseq = models.CharField(null = True, blank = True, max_length = 3000)
    fail = models.CharField(null = True, blank = True, max_length = 500)
    failures = models.IntegerField(default = 0, blank = True)
    user = models.ForeignKey(User, null = True, blank = True, on_delete = models.CASCADE)
    def __unicode__(self):
        return self.name