 # -*- coding: utf-8 -*-
"""
@author: aitana neves | heig-vd
"""

global find_all
global smooth
global smooth_wo_boundary
global get_continuous_fragments
global list2string_wo_spaces
global indices2string
global protter_plot_struct
global protter_plot_conservation
global new_figure
global jaccard_shoulder

#%% IMPORTS

import numpy as np
import urllib2  
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#%% GLOBAL FUNCTIONS

def find_all(a_str, sub):
    
    start = 0
    idx= []
    while start >= 0:
        start = a_str.find(sub, start)
        idx.append(start)
        if start == -1:
            break
        start += len(sub)
    return idx[0:-1]

def smooth(y, box_pts):
    
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
    
def smooth_wo_boundary(y, box_pts):

    box = np.ones(box_pts)/box_pts
    y_smooth = (np.convolve(y, box, mode='same')).tolist()
    r = int(round(box_pts))
    yf = []
    yf = list(y[0:r])
    yf.extend(y_smooth[r:-r+1])
    yf.extend(y[-r:-1])
    return yf     

def get_continuous_fragments(input):
    
    indices2 = (np.where(np.diff(input)>1)[0]).tolist()
    continuous_fragments = []
    indices2[:0] = [-1] # prepend a -1
    indices2.extend([len(input)-1]) # append len()
    reduce(lambda x, y: continuous_fragments.append(input[x+1:y+1]) or y, indices2)
    continuous_fragments = filter(None,continuous_fragments) # remove empty lists (consisting of only one amino acid) 
    return continuous_fragments
    
def protter_plot_struct(path,short_name,protein_name,version,indices_A,indices_E,indices_DR,indices_B):
      
    protter_query = 'http://wlab.ethz.ch/protter/create?up='+ protein_name + \
                        '&nterm=intra&tm=Phobius.TM,UP.TRANSMEM' + \
                        '&mc=lightgray&lc=blue&tml=numcount&numbers&legend' + \
                        '&n:SIGNAL-PEPTIDE,cc:blue=UP.SIGNAL' + \
                        '&n:ALPHA-HELIX,bc:red=' + indices_A + \
                        '&n:BETA-STRANDS,bc:chartreuse=UP.STRAND,' + indices_E + \
                        '&n:DISORDERED,fc:cornflowerblue=' + indices_DR + \
                        '&n:BURIED,s:box,fc:mediumorchid=' + indices_B + \
                        '&n:PTMs,s:diamond,fc:orange=UP.CARBOHYD,UP.MOD_RES,UP.CROSSLNK' + '&format=pdf'
    response = urllib2.urlopen(protter_query)
    protter_file = open(path+short_name+'_PROTTER-STRUCTURE'+version+'.pdf','w')
    protter_file.write(response.read())
    protter_file.close()
        
def list2string_wo_spaces(mylist):
    
    output_string = (((str(mylist)).replace(', ',',')).replace("'",'')).strip('[]')
    return output_string
            
def indices2string(indices):

    indices_string = (str((np.array(indices) + 1).tolist()).strip('[]')).replace(' ','')
    return indices_string
            
def protter_plot_conservation(path,short_name,protein_name,version,indices_VAR,indices_VAR_BLOSUM,indices_CONS_NR,indices_SC1,indices_SC):
              
    protter_query2 = 'http://wlab.ethz.ch/protter/create?up='+protein_name + \
                        '&nterm=intra&tm=Phobius.TM,UP.TRANSMEM' + \
                        '&mc=lightgray&lc=blue&tml=numcount&numbers&legend' + \
                        '&n:ISOLATES-VARIANTS,s:circ,fc:magenta=UP.VARIANT,' + indices_VAR + \
                        '&n:ISOLATES-VARIANTS-BLOSUM95,s:circ,cc:magenta=' + indices_VAR_BLOSUM + \
                        '&n:NR-CONSERVED-BLOSUM95,s:box,fc:green=' + indices_CONS_NR + \
                        '&n:EVOL_VAR,s:circ,bc:skyblue=' + indices_SC1 + \
                        '&n:EVOL-CONS,s:circ,bc:firebrick=' + indices_SC + \
                        '&n:REPEATS,s:diamond=UP.REPEAT' + '&format=pdf'   
    response2 = urllib2.urlopen(protter_query2)
    protter_file2 = open(path+short_name+'_PROTTER-CONSERVATION.pdf','w')
    protter_file2.write(response2.read())
    protter_file2.close()
    
def new_figure(size=(4,3)):
    
    fig = plt.figure(num=None,figsize=size,dpi=72)
    ax = fig.add_axes([0.17,0.02,0.7,0.62])
    return(fig,ax)
    
def jaccard_shoulder(list1,list2,shoulder=0.5):
    
    if type(list1) == np.ndarray:
        list1 = list1.tolist()
    if type(list2) == np.ndarray:
        list2 = list2.tolist()
    
    # add a "shoulder" to boolean list1 and list2 (i.e. add 0.5 to a 0-valued pixel neighboring a 1-valued pixel)

    max_envelope = []
    min_envelope = []
    list1s = list(list1)    
    if (list1[1])==1 & (list1[0]==0):
        list1s[0] = shoulder
    if (list1[-2]==1) & (list1[-1]==0):
        list1s[-1] = shoulder
    
    list2s = list(list2)    
    if (list2[1]==1) & (list2[0]==0):
        list2s[0] = shoulder
    if (list2[-2]==1) & (list2[-1]==0):
        list2s[-1] = shoulder
    
    max_envelope.append(max(list1s[0],list2s[0]))  
    min_envelope.append(min(list1s[0],list2s[0])) 
    
    for i in range(1,len(list1)-1):    
        if ((list1[i+1]==1) | (list1[i-1]==1)) & (list1[i]==0):
            list1s[i] = shoulder
        if ((list2[i+1]==1) | (list2[i-1]==1)) & (list2[i]==0):
            list2s[i] = shoulder
        max_envelope.append(max(list1s[i],list2s[i]))  
        min_envelope.append(min(list1s[i],list2s[i])) 
    
    max_envelope.append(max(list1s[-1],list2s[-1]))  
    min_envelope.append(min(list1s[-1],list2s[-1])) 
    
    # compute Jaccard-like distance using list1s and list2s (i.e. min-envelope/max-envelope)   
    max_envelope = sum(max_envelope)
    min_envelope = sum(min_envelope)
    if max_envelope==0:
        return 0
    d = min_envelope/float(max_envelope)    
    
    return d
