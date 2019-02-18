# -*- coding: utf-8 -*-
"""
@author: aitana neves | heig-vd
"""

global get_uniprot_info
global get_disordered_regions_PrDOS
global get_secondary_structure_JPred
global get_scratch_predictions
global get_solvent_accessibility_JPred
global get_hydrophobicity
global get_consensus_predictions
global write_protein_characterization_to_file

#%% IMPORTS

import pandas as pd    
import numpy as np

#%% GLOBAL FUNCTIONS

def get_uniprot_info(uniprot_ref):
    
    print('\n>> PROTEIN CHARACTERIZATION | RETRIEVING UNIPROT INFO\n')
    
    from xml.dom import minidom
    import urllib
    
    url = 'http://www.uniprot.org/uniprot/'+uniprot_ref+'.xml'
    dom = minidom.parse(urllib.urlopen(url))
    PTM = [] # post-translational modifications (glycolysation sites)
    mutations = []
    sTM = []
    eTM = []
    sSP = []
    eSP = []        
    for node in dom.getElementsByTagName('feature'):
        type = str(node.attributes['type'].value)
        if type=='transmembrane region':
            sTM = node.getElementsByTagName('begin')
            sTM = int(sTM[0].attributes['position'].value)
            eTM = node.getElementsByTagName('end')
            eTM = int(eTM[0].attributes['position'].value)
            
        if type=='signal peptide':
            sSP = node.getElementsByTagName('begin')
            sSP = int(sSP[0].attributes['position'].value)
            eSP = node.getElementsByTagName('end')
            eSP = int(eSP[0].attributes['position'].value)
            
        if type=='glycosylation site':
            ptm = node.getElementsByTagName('position')
            PTM.append(int(ptm[0].attributes['position'].value))
        if type =='mutagenesis site':
            mut = node.getElementsByTagName('position')
            mutations.append(int(mut[0].attributes['position'].value))
    try:
        indices_TM = range(sTM-1,eTM,1)       # ! uniprot starts counting from 1, not 0 !
    except:
        indices_TM = []
    try:
        indices_SIGNAL = range(sSP-1,eSP,1)   # ! uniprot starts counting from 1, not 0 !
    except:
        indices_SIGNAL = []            
    indices_PTM = (np.array(PTM)-1).tolist()  # ! uniprot starts counting from 1, not 0 !
    indices_mut = (np.array(mutations)-1).tolist()
    return (indices_TM,indices_SIGNAL,indices_PTM,indices_mut)

def get_disordered_regions_PrDOS(data_folder):

    dr = pd.read_csv(data_folder+'PrDOS.csv',sep=',',header=None)
    indices_DR = np.where(dr[2]==1)
    indices_DR = np.array(indices_DR)
    indices_DR = indices_DR + 1
    indices_DR= (str(indices_DR.tolist())).strip('[]')
    indices_DR = indices_DR.replace(' ','')
    
    return indices_DR
    
def get_secondary_structure_JPred(data_folder):
    
    fo = open(data_folder+'JPred.txt','r')   
    dd = fo.readlines()
    secondary_structure = dd[1][9:]
        
    secondary_structure = secondary_structure.replace('-','0')
    secondary_structure = secondary_structure.replace('H','1')
    secondary_structure = secondary_structure.replace('E','2')
    secondary_structure = secondary_structure[0:-2]
    secondary_structure = secondary_structure.split(',')
    secondary_structure = np.array(secondary_structure)
    indices_H = np.where(secondary_structure=='1')
    indices_H = np.array(indices_H)
    indices_H = indices_H + 1
    indices_H = (str(indices_H.tolist())).strip('[]')
    indices_H = indices_H.replace(' ','')
    indices_E = np.where(secondary_structure=='2')     
    indices_E = np.array(indices_E)
    indices_E = indices_E + 1
    indices_E = (str(indices_E.tolist())).strip('[]')
    indices_E = indices_E.replace(' ','')  

    return (indices_H,indices_E)     
    
def get_solvent_accessibility_JPred(data_folder):

    fo = open(data_folder+'JPred.txt','r')   
    dd = fo.readlines()
    
    buried = dd[4][9:];    
    buried = buried.replace('-','0')
    buried = buried.replace('B','5')
    buried = buried[0:-2]
    buried = buried.split(',')
    buried = np.array(buried)

    buried0 = dd[5][9:].replace('-','0').replace('B','1').split(',')[0:-1]      
    buried5 = dd[4][9:].replace('-','0').replace('B','1').split(',')[0:-1]  
    buried25 = dd[3][10:].replace('-','0').replace('B','1').split(',')[0:-1]  
    
    buried_score = np.array([float(buried0[i]) for i in range(len(buried0))]) + np.array([float(buried5[i]) for i in range(len(buried5))]) + np.array([float(buried25[i]) for i in range(len(buried25))])
    buried_score[np.where(buried_score==0)] = float('nan')
    buried_score[np.where(buried_score==3)] = 0
    buried_score[np.where(buried_score==2)] = 3 #5
    buried_score[np.where(buried_score==1)] = 2 #25  
    buried_score[np.where(buried_score==3)] = 1 #only way to rename numbers in buried_score (need an "extra container")     
    
    indices_B = np.where(buried=='1')
    indices_B = np.array(indices_B)
    indices_B = indices_B + 1
    indices_B = ((str(indices_B.tolist())).strip('[]')).replace(' ','') 
    
    return (indices_B,buried_score)
    
def get_scratch_predictions(data_folder):

    fo = open(data_folder+'SCRATCH.txt','r')     
    fo.seek(0)
    dd = fo.readlines()   
    ss = np.array(list((dd[1]).strip('\n')))
    indices_H_scratch = ((str((np.array(np.where(ss=='H')) + 1).tolist())).strip('[]')).replace(' ','')
    indices_E_scratch = ((str((np.array(np.where(ss=='E')) + 1).tolist())).strip('[]')).replace(' ','')
    disordered_scratch = np.array(list((dd[2]).strip('\n')))
    indices_DR_scratch = ((str((np.array(np.where(disordered_scratch=='D')) + 1).tolist())).strip('[]')).replace(' ','')
    buried_scratch = np.array(list((dd[3]).strip('\n')))
    indices_B_scratch = ((str((np.array(np.where(buried_scratch=='-')) + 1).tolist())).strip('[]')).replace(' ','')
    
    return (indices_H_scratch,indices_E_scratch,indices_DR_scratch,indices_B_scratch)
    
def get_hydrophobicity(Hphob_scale,refseq):

    position_specific_hydrophilicity_score = np.array([-float(Hphob_scale[refseq[i]]) for i in range(len(refseq))])
    
    return position_specific_hydrophilicity_score

def get_consensus_predictions(list1,list2):
    
    indices_consensus = list(set(list1.split(',')) & set(list2.split(',')))
    
    return indices_consensus
    
def write_protein_characterization_to_file(folder,L,indices_H_consensus,indices_E_consensus,indices_PTM):
    
    fo = open(folder+'Position-specific_secondary-structure_PTM.txt','w')
    for i in range(L):
        tmp = str(i+1)+'\t'
        if indices_H_consensus.count(str(i+1)) > 0:
            tmp += 'alpha\t'
        elif indices_E_consensus.count(str(i+1)) > 0:
            tmp += 'beta\t'
        else:
            tmp += 'coil\t'
        if indices_PTM.count(i) > 0:
            tmp += 'PTM'
        tmp += '\n'
        fo.write(tmp)
    fo.close()
    