# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 16:14:33 2018

@author: jibril
"""

from django import template

register = template.Library()

@register.filter
def modulo(num, val):
    return num % val
