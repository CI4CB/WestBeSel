# -*- coding: utf-8 -*-
"""
@author: aitana
"""

global get_current_directory
global set_directories

#%% IMPORTS

import os
from epitope_prediction import useful_functions as uf

#%% GLOBAL FUNCTIONS

def get_current_directory(main_file):

    main_path = os.path.realpath(main_file)
    slash_idx = uf.find_all(main_path,'/')
    current_directory = main_path[0:slash_idx[-1]+1]

    return current_directory

def set_directories(wdir,project_name,sf):
    
    import os
    import collections
    
    print('\n>> DIRECTORIES | Creating folders for '+project_name+' (if not exist):\n')    
    
    if wdir[-1]!='/': # in case the user forgot the slash at the end
        wdir = wdir+'/'
    mypathR = wdir+project_name+r'/'
    if not os.path.isdir(mypathR):
        os.mkdir(mypathR)
    print(mypathR+'\n')

    subf = collections.defaultdict(dict)
    for i in range(len(sf)):
        tmp_path = mypathR+sf[i]+'/'
        subf[sf[i]] = {'path':tmp_path,'index':i}        
        if not os.path.isdir(tmp_path):
            os.mkdir(tmp_path)
        print('..... '+str(i)+'_'+sf[i]+'\n')
        
    return (mypathR,subf)