# -*- coding: utf-8 -*-
"""
@author: aitana neves | heig-vd
"""

global get_all_epitope_predictions
global view_epitope_predictions
global get_consensus_epitopes
global crop_epitopes
global post_process_cysteins
global terminal_epitopes
global discard_epitopes_with_repeats
global write_epitope_ws_fasta
global draw_epitope_ws_weblogos
global blast_all_epitopes
global get_all_sliding_windows
global keep_epitopes_with_conserved_hydrophilic_stretches
global compute_pareto_dominance
global discard_epitopes_with_PTM
global discard_hydrophobic_epitopes

#%% IMPORTS

import os
import epitope_prediction as ep
from epitope_prediction import useful_functions as uf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import collections
import pickle
import pandas as pd
from matplotlib.ticker import MultipleLocator
import matplotlib
import scipy.cluster.hierarchy as cl
import scipy.spatial.distance as ssd
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SearchIO
import xlsxwriter
import pickle

#%% LOCAL FUNCTIONS

def visualize_epitopes(ml,score_res,folder,project_name,met):

    plt.rcParams['grid.linewidth'] = 0.2  
    fig = plt.figure(num=None,figsize=(5,3),dpi=72)
    ax = fig.add_axes([0.15,0.15,0.75,0.75])
    plt.plot(score_res,'b-')
    plt.xlabel(r'# amino acid')
    plt.ylabel(r'Epitopes')
    plt.title(project_name)
    ax.xaxis.grid() # vertical lines       
    ax.xaxis.set_minor_locator(MultipleLocator(25))
    ax.xaxis.grid(True,'minor',linewidth=0.5)
    plt.ylim(0, 1.5)
    plt.xlim(0,len(score_res))
    plt.hlines(y=0, xmin=0, xmax=len(score_res),linestyle='dotted')
    yloc = [1.28,1.3,1.32]
    for kdx in range(len(ml)):
        loc = yloc[kdx%len(yloc)]
        plt.hlines(y=loc,xmin=ml[kdx][0],xmax=ml[kdx][-1],linewidth=1.5,linestyle='solid')   
    #plt.show()
    fig.savefig(folder+met+'_visualization.pdf', dpi=fig.dpi)

def AAP_epitopes(data_folder,L):

    print '..... AAP'
    fo = open(data_folder+'AAP_20aa.txt','r') 
    fo.seek(0)
    d0 = fo.readlines()
    if len(d0)>0:
        # get most likely and least likely epitope fragments (peptides, 5-18 residues)
        df_20 = pd.read_csv(data_folder+'AAP_20aa.txt',sep='\t',header=None)
        df_20_NO = pd.read_csv(data_folder+'AAP_20aa_NO.txt',sep='\t',header=None)
        df_12_NO = pd.read_csv(data_folder+'AAP_12aa_NO.txt',sep='\t',header=None)
        df_16_NO = pd.read_csv(data_folder+'AAP_16aa_NO.txt',sep='\t',header=None)
        L20 = 20
        L12 = 12
        L16 = 16
    
        # determine residue positions for each fragment (overlapping)
        ml1 = [] 
        positives = []
        positive_score = []
        for kdx in range(len(df_20)):#[1:-1]:
            ml1.append(range(df_20[0][kdx],df_20[0][kdx]+L20))
            positives.extend(range(df_20[0][kdx],df_20[0][kdx]+L20))      
            positive_score.extend((np.ones(L20,dtype=np.float)*df_20[2][kdx]).tolist())
        positive_score = np.array(positive_score)
        
        # for each residue, calculate score
        score_res = []        
        for kdx in range(L):
            indices_pos = [i for i, x in enumerate(positives) if x == kdx] # index of positive epitope fragments where the residue appears
            if (kdx < L20-1):   # !! correct for "boundary effects" due to the 20aa sliding window
                number_of_possible_positives = kdx+1
            elif (kdx > L-L20-1):
                number_of_possible_positives = L-kdx        
            else:
                number_of_possible_positives = L20
            score_res.append(sum(positive_score[indices_pos])/float(number_of_possible_positives)) # positive scores of the epitope fragments where the residue appears
    
        score_res = np.array(score_res)   
        score_res = score_res/max(score_res)    
    
        # non overlapping
        ml_NO = []
        for kdx in range(len(df_20_NO)):
            ml_NO.append(range(df_20_NO[0][kdx],df_20_NO[0][kdx]+L20))
            
        ml12_NO = []
        for kdx in range(len(df_12_NO)):
            ml12_NO.append(range(df_12_NO[0][kdx],df_12_NO[0][kdx]+L12))
            
        ml16_NO = []
        for kdx in range(len(df_16_NO)):
            ml16_NO.append(range(df_16_NO[0][kdx],df_16_NO[0][kdx]+L16))
    else:
        return([],[])
        # keep regions with a positive score only
    #input = (np.where(score_res>0)[0]).tolist()
    #ml = uf.get_continuous_fragments(input)
    #ml_O = ml
    
    ml = ml_NO
    
    return (ml,score_res)
    
def ABCPred_epitopes(data_folder,L):

    print '..... ABCPred'

    df = pd.read_csv(data_folder+'ABCPred.txt',sep='\t',header=None)   
    window = len(df[1][0])
    ml = []  
    score_res = []
    positives = []
    positive_score = []
    perc = 0.2 # 20% best
    for kdx in range(int(round(perc*len(df)))): # keep best 20%
        ml.append(range(df[2][kdx],df[2][kdx]+window))
        positives.extend(range(df[2][kdx],df[2][kdx]+window))
        positive_score.extend((np.ones((window,1),dtype=np.float)*df[3][kdx]).tolist())  
    positive_score = np.array(positive_score)
        
    score_res = []
    for kdx in range(L):
        try:
            indices_pos = [i for i, x in enumerate(positives) if x == kdx] # index of positive epitope fragments where the residue appears                         
            score_res.append(sum(positive_score[indices_pos])) # positive scores of the epitope fragments where the residue appears
        except:
            score_res.append(0)
    score_res = np.array(score_res)
    score_res = (score_res-min(score_res))/max(score_res-min(score_res))    
        
    score_res = uf.smooth(score_res.tolist(),5)
    
    # keep residues with score_res > 0.5
    perc2 = 0
    input = (np.where(score_res>perc2)[0]).tolist()
    ml = uf.get_continuous_fragments(input)

    return (ml,score_res)
        
def BCPred_epitopes(data_folder,L):
    
    print '..... BCPred'    
    fo = open(data_folder+'BCPred_20aa.txt','r') 
    fo.seek(0)
    d0 = fo.readlines()
    if len(d0)>0:
        df_20 = pd.read_csv(data_folder+'BCPred_20aa.txt',sep='\t',header=None)
        df_20_NO = pd.read_csv(data_folder+'BCPred_20aa_NO.txt',sep='\t',header=None)
        df_12_NO = pd.read_csv(data_folder+'BCPred_12aa_NO.txt',sep='\t',header=None)
        df_16_NO = pd.read_csv(data_folder+'BCPred_16aa_NO.txt',sep='\t',header=None)
        L20 = 20
        L12 = 12
        L16 = 16
    
        # determine residue positions for each fragment
        ml1 = [] 
        positives = []
        positive_score = []
        for kdx in range(len(df_20)):#[1:-1]:
            ml1.append(range(df_20[0][kdx],df_20[0][kdx]+L20))
            positives.extend(range(df_20[0][kdx],df_20[0][kdx]+L20))      
            positive_score.extend((np.ones(L20,dtype=np.float)*df_20[2][kdx]).tolist())
        positive_score = np.array(positive_score)
        
        # for each residue, calculate score
        score_res = []        
        for kdx in range(L):
            indices_pos = [i for i, x in enumerate(positives) if x == kdx] # index of positive epitope fragments where the residue appears
            if (kdx < L20-1):   # !! correct for "boundary effects" due to the 20aa sliding window
                number_of_possible_positives = kdx+1
            elif (kdx > L-L20-1):
                number_of_possible_positives = L-kdx        
            else:
                number_of_possible_positives = L20
            score_res.append(sum(positive_score[indices_pos])/float(number_of_possible_positives)) # positive scores of the epitope fragments where the residue appears
    
        score_res = np.array(score_res)   
        score_res = score_res/max(score_res)    
    
        # non overlapping
        ml_NO = []
        for kdx in range(len(df_20_NO)):
            ml_NO.append(range(df_20_NO[0][kdx],df_20_NO[0][kdx]+L20))
            
        ml12_NO = []
        for kdx in range(len(df_12_NO)):
            ml12_NO.append(range(df_12_NO[0][kdx],df_12_NO[0][kdx]+L12))
            
        ml16_NO = []
        for kdx in range(len(df_16_NO)):
            ml16_NO.append(range(df_16_NO[0][kdx],df_16_NO[0][kdx]+L16))
    else:
        return([],[])
    # keep regions with a positive score only
    #input = (np.where(score_res>0)[0]).tolist()
    #ml = a.get_continuous_fragments(input)
    #ml_O = ml
    
    ml = ml_NO        
    
    return (ml,score_res)
    
def BepiPred_epitopes_one(data_folder,L):
    
    print '..... BepiPred'    
    
    # epitopic residues
    fo = open(data_folder+'BepiPred.txt','r') 
    fo.seek(0)
    d0 = fo.readlines()
    if len(d0)>0:
        score_res = [float(d0[i].split(' ')[-7]) for i in range(len(d0))]
        score_res = np.array(score_res)
        score_res[score_res>0] = score_res[score_res>0]/max(score_res)
        score_res[score_res<0] = -score_res[score_res<0]/min(score_res)
    
        input = [i for i in range(len(d0)) if d0[i].split('|')[1].strip('\n')=='E'] # find positions where there is an "E" instead of '.', i.e. epitopes
        ml = uf.get_continuous_fragments(input)
        
        return (ml,score_res)
    else:
        return([],[])

def BepiPred_epitopes(data_folder,L):
    
    print '..... BepiPred 2.0'    
    
    # epitopic residues
    fo = open(data_folder+'BepiPred.txt','r') 
    fo.seek(0)
    d0 = fo.readlines()
    if len(d0)>0:
        score_res = [float(d0[i].split('\t')[-1]) for i in range(1,len(d0))]
        score_res = np.array(score_res)
        ml = [i if score_res[i]>0.5 else 0  for i in range(len(score_res))]
        ml = uf.get_continuous_fragments(ml)
        return (ml,score_res)
    else:
        return([],[])
    
def CBTope_epitopes(data_folder,L):
    
    print '..... CBTope'    
    
    THRESHOLD = ep.CBTope_THRESHOLD
    
    # scores (above 4 considered as epitope)
    fo = open(data_folder+'CBTope_scores.txt','r')   
    d0 = fo.read()
    d1 = (d0.replace('&nbsp; ','')).strip('\n')
    score_res0 = list(d1)   
    score_res = np.array(score_res0,dtype=float)/4.5-0.5  

    ml = np.array(score_res0,dtype=int)
    ml[ml<THRESHOLD] = 0
    ml[ml>0] = 1
    input = (np.where(ml==1)[0]).tolist()
    ml = uf.get_continuous_fragments(input)
    
    return (ml,score_res)
    
def COBEPro_epitopes(data_folder,L):
    
    print '..... COBEPro'    
    
    # get most likely and least likely epitope fragments (peptides, 5-18 residues)
    ml0 = pd.read_csv(data_folder+'COBEPro_ML.txt',sep=' ',header=None)
    ll0 = pd.read_csv(data_folder+'COBEPro_LL.txt',sep=' ',header=None)
    
    mlL = [] # length of epitope fragments (most likely)
    llL = [] # length of epitope fragments (least likely)
    for i in range(len(ml0[2])): mlL.append(len(ml0[2][i]))     
    for i in range(len(ll0[2])): llL.append(len(ll0[2][i]))   
    
    # determine residue propensity score based on epitope fragments score (see Methods COBEPro paper)    
    # identify top 5% and bottom 5% from ml and ll, respectively
    perc = 1#0.75
    top = range(int(round(perc*len(ml0))))
    bottom = range(int(round(perc*len(ll0))))
    ml0[0] = ml0[0][top]
    ml0[2] = ml0[2][top]
    ml0[1] = ml0[1][top]
    ll0[0] = ll0[0][bottom]
    ll0[2] = ll0[2][bottom]
    ll0[1] = ll0[1][bottom]
    
    # determine residue positions for each fragment
    ml1 = []        
    positives = []
    positive_score = []
    for kdx in range(len(top)):
        tmp = range(int(ml0[1][kdx]),int(ml0[1][kdx]+len(ml0[2][kdx])))
        ml1.append(tmp)
        positives.extend(tmp)        
        positive_score.extend((np.ones((len(ml0[2][kdx]),1),dtype=np.float)*ml0[0][kdx]).tolist())
    positive_score = np.array(positive_score)
    
    ll1 = []        
    negatives = []
    negative_score = []
    for kdx in range(len(bottom)):
        tmp = range(int(ll0[1][kdx]),int(ll0[1][kdx]+len(ll0[2][kdx])))
        ll1.append(tmp)
        negatives.extend(tmp) 
        negative_score.extend((np.ones((len(ll0[2][kdx]),1),dtype=np.float)*ll0[0][kdx]).tolist())
    negative_score = np.array(negative_score) 
    
    # for each residue, calculate score
    score_res = []        
    for kdx in range(L):
        indices_pos = [i for i, x in enumerate(positives) if x == kdx] # index of positive epitope fragments where the residue appears
        indices_neg = [i for i, x in enumerate(negatives) if x == kdx] # index of negative epitope fragments where the residue appears
        score_res.append(sum(positive_score[indices_pos]) + sum(negative_score[indices_neg])) # positive scores of the epitope fragments where the residue appears

    score_res = np.array(score_res)   
    score_res[score_res>0] = score_res[score_res>0]/max(score_res)
    score_res[score_res<=0] = -score_res[score_res<=0]/min(score_res)
    
    # keep regions with a positive score only
    input = (np.where(score_res>0)[0]).tolist()
    ml = uf.get_continuous_fragments(input)
    
    return (ml,score_res)
    
def SVMTriP_epitopes(data_folder,L):

    print '..... SVMTriP'    
    
    score_res = np.zeros((1,L),dtype=float)[0].tolist()
    try: 
        # 20aa
        df = pd.read_csv(data_folder+'SVMTriP.txt',sep='\t',header=None)   
        ml = []        
        for kdx in range(len(df)):
            ml.append(range(df[1][kdx],df[2][kdx]+1))
    except:
        ml=[]
        score_res=np.array(score_res)

    return (ml,score_res)
    
def TEPRF_epitopes(data_folder,L):
    
    print '..... TEPRF'        
    
    df = pd.read_csv(data_folder+'TEPRF.txt',sep='\t',header=None)   
    
    if L != len(df) + 20 - 1:
        print '..... !!!!! PREDICT_EPITOPES.TEPRF_epitopes: something weird !!!!!'
        print df
     
    # identify top 20% and bottom 20% from ml and ll, respectively
    perc = 0.2
    ind_pos = np.where(np.array(df[3])>0) 
    ml0 = (df[4][ind_pos[0].tolist()]).tolist() # positive scores     
    ind_neg = np.where(np.array(df[3])<0) 
    ll0 = (df[4][ind_neg[0].tolist()]).tolist() # negative scores
    ids_pos = np.argsort(np.array(ml0)) # sort scores (! ascending order!)
    ind_pos = ind_pos[0][ids_pos]
    ids_neg = np.argsort(np.array(ll0)) # sort scores
    ind_neg = ind_neg[0][ids_neg]
    top = ind_pos[int(round((1-perc)*len(ind_pos))):len(ind_pos)]
    bottom = ind_neg[int(round((1-perc)*len(ind_neg))):len(ind_neg)]
    
    ml1 = []        
    positives = []
    positive_score = []
    for kdx in top:
        tmp = range(kdx,kdx+20)
        ml1.append(tmp)
        positives.extend(tmp)        
        positive_score.extend((np.ones(20,dtype=np.float)*df[4][kdx]).tolist())
    positive_score = np.array(positive_score)
    
    ll1 = []        
    negatives = []
    negative_score = []
    for kdx in bottom:
        tmp = range(kdx,kdx+20)
        ll1.append(tmp)
        negatives.extend(tmp) 
        negative_score.extend((np.ones(20,dtype=np.float)*df[4][kdx]).tolist())
    negative_score = np.array(negative_score) 
    
    # for each residue, calculate score
    score_res = []        
    for kkdx in range(L):
        indices_pos = [i for i, x in enumerate(positives) if x == kkdx] # index of positive epitope fragments where the residue appears
        indices_neg = [i for i, x in enumerate(negatives) if x == kkdx] # index of negative epitope fragments where the residue appears
        score_res.append(sum(positive_score[indices_pos]) - sum(negative_score[indices_neg])) # positive scores of the epitope fragments where the residue appears

    score_res = np.array(score_res)   
    score_res[score_res>0] = score_res[score_res>0]/max(score_res)
    score_res[score_res<=0] = -score_res[score_res<=0]/min(score_res)               
    
    # keep regions with a positive score only
    input = (np.where(uf.smooth(score_res,5)>0)[0]).tolist()
    ml = uf.get_continuous_fragments(input)
    
    return (ml,score_res)
    
def Epitopia_epitopes(data_folder,L):
    
    print '..... Epitopia'    
    
    ## epitopic residues
    df = pd.read_csv(data_folder+'Epitopia.txt',sep='\t',header=None)
    indices = [int(df[0][i][3:]) for i in range(1,len(df))] # identify indices to which scores are associated
    sid = np.argsort(np.array(indices))
    indices = np.array(indices)[sid]   
    score_res = np.array(df[2][1:],dtype=float)
    score_res = score_res[sid]
    score_res = (score_res-min(score_res))/max(score_res-min(score_res))      
    score_res = uf.smooth(score_res.tolist(),3)

    #perc = 0. #0.25
    #input = np.where(score_res>perc)[0].tolist()
    #ml = uf.get_continuous_fragments(input)

    ## best epitopes as predicted by epitopia
    fo = open(data_folder+'Epitopia_best-epitopes.txt','r') 
    fo.seek(0)
    d0 = fo.readlines()
    ind_plus = [i for i in range(len(d0)) if d0[i]=='++++++++++++++++++++++++++\n'] # find lines with '+'
    ind_plus[:0] = [-1]        
    st = [int(d0[ind_plus[i]+5].strip('\n')[3:])-1 for i in range(len(ind_plus)-1)]
    ea = [int(d0[ind_plus[i+1]-2].strip('\n')[3:])-1 for i in range(len(ind_plus)-1)]
    epitopia_ml = [range(st[i],ea[i]+1) for i in range(len(st))]
    
    ## identify 8 best epitopes by averaging residue scores
    #av_score_ml = [sum(score_res[ml[i]])/len(ml[i]) for i in range(len(ml))] 
    #sid2 = np.argsort(av_score_ml)
    #best_ml = []
    #best_ml = [ml[sid2[i]] for i in range(-9,-1)]
    
    #ml_all = ml
    ml = epitopia_ml   
    
    return (ml,score_res)
    
def plot_matrix(folder,c2,labels):
    
    N_methods = len(ep.methods)
    fig = plt.figure(num=None,figsize=(4,3),dpi=300)
    ax = fig.add_axes([0.17,0.02,0.7,0.62])
    ax.xaxis.set_tick_params(labeltop='on',labelbottom='off')
    plt.xticks(range(0,N_methods))
    ax.set_xticklabels(labels,rotation=90)
    plt.yticks(range(0,N_methods))
    ax.set_yticklabels(ep.methods,rotation=0)
    plt.imshow(c2,interpolation='nearest',cmap=plt.cm.coolwarm)
    plt.colorbar()
#    plt.show()
    plt.savefig(folder+'FIGURE-hierarchical_clustering_methods.pdf', dpi=fig.dpi)
    
def remove_diag(mat,val=0):
    
    mat2 = np.array(mat)
    for i in range(np.shape(mat)[0]):
        mat2[i,i] = val
    return mat2
    
def check_overlap(ml,indices_SIGNAL,indices_TM,indices_PTM):
    
    nb_epitopes = len(ml)
    overlap_SP = [ list(set(ml[i]).intersection(indices_SIGNAL)) for i in range(nb_epitopes) ]
    overlap_TM = [ list(set(ml[i]).intersection(indices_TM)) for i in range(nb_epitopes) ]
    overlap_PTM = [ list(set(ml[i]).intersection(indices_PTM)) for i in range(nb_epitopes) ]
    mec= []
    for kdx in range(nb_epitopes): # a fragment could overlap many protein domains... still, only last will be considered for colouring
        if len(overlap_SP[kdx])>0:
            mec.append([0,1,0,1])
        elif len(overlap_TM[kdx])>0:
            mec.append([1,0.55,0,1])
        elif len(overlap_PTM[kdx])>0:
            mec.append([1,0,1,1])
        else:
            mec.append([0,0,0,1])        
    return mec
    
def compute_pareto_dominance(data):
    try:
        nb_features = np.shape(data)[1]   
        nb_epitopes = np.shape(data)[0]      
        pareto_dominance = np.zeros((nb_epitopes,nb_epitopes),dtype=int)
    except:
        return (),()
    for i in range(nb_epitopes):
        for j in range(nb_epitopes): #range(i+1,nb_epitopes,1):
            greater_or_equal = [k for k in range(nb_features) if data[i][k]>=data[j][k]] # check each feature
            greater = [k for k in range(nb_features) if data[i][k]>data[j][k]]
            pareto_dominance[i][j] = 1 if ( (len(greater_or_equal)==nb_features) & (len(greater)>=1) ) else 0
    pareto_dominant = np.sum(pareto_dominance,1).tolist()
    idx_pareto = np.argsort(pareto_dominant)[::-1] # [::-1] is to take the reverse order, to have sorting in descending values rather than ascending              
    pareto_dominant = np.array(pareto_dominant)/float(nb_epitopes-1) # normalize Pareto score by nb_epitopes 
    
    return (idx_pareto,pareto_dominant)
    
def write_excel_file(folder,project_name,ml,idx_pareto,mec,refseq,av_score_isolates,hydrophilicity_score,av_score,av_solvent_accessibility,pareto_dominant,max_repeat_length):
    
    workbook = xlsxwriter.Workbook(folder+'LIST-final-epitopes.xls')
    format_heading = workbook.add_format({'bold': True})
    format_bold = workbook.add_format({'bold': True})
    format_cell_color = workbook.add_format()
    format_bold.set_text_wrap()
    format_cell = workbook.add_format()
    format_bold.set_font_size(9)
    format_cell.set_font_size(7)
    format_cell.set_align('center')
    format_cell.set_align('vcenter')       
    format_bold.set_align('center')
    format_bold.set_align('vcenter')  
    format_cell_color.set_font_size(7)
    format_cell_color.set_align('center')
    format_cell_color.set_align('vcenter')       
    format_cell_color.set_font_color('magenta')
    format_cell_color.set_italic()
    worksheet = workbook.add_worksheet()
    worksheet.set_landscape()
    worksheet.set_column('A:A', 4)
    worksheet.set_column('B:B', 3)
    worksheet.set_column('C:C', 10)
    worksheet.set_column('D:J', 5)
    worksheet.merge_range('A1:C1',project_name,format_heading)
    worksheet.write(1,0,'Start',format_bold)
    worksheet.write(1,1,'End',format_bold)
    worksheet.write(1,2,'Fragment',format_bold)
    worksheet.write(1,3,'Within-isolate\nconservation',format_bold)
    worksheet.write(1,4,'Average\nhydrophilicity',format_bold)
    worksheet.write(1,5,'Solvent accessibility',format_bold)
    worksheet.write(1,6,'jury-score',format_bold)
    worksheet.write(1,7,'Pareto-score',format_bold)
    worksheet.write(1,8,'Repeats',format_bold)
    worksheet.write(1,9,'xml Number',format_bold)
    for x,i in enumerate(idx_pareto):
        worksheet.write(x+2,0,ml[i][0]+1,format_cell)  # add +1 to start counting from 1 (more intuite in xls file and similar to Blast output)
        worksheet.write(x+2,1,ml[i][-1]+1,format_cell) # add +1 to start counting from 1 (more intuite in xls file and similar to Blast output)
        if mec[i]==[1,0,1,1]:
            worksheet.write(x+2,2,refseq[ml[i][0]:ml[i][-1]+1],format_cell_color)
        else:
            worksheet.write(x+2,2,refseq[ml[i][0]:ml[i][-1]+1],format_cell)
        try:
            worksheet.write(x+2,3,round(av_score_isolates[i],2),format_cell)
        except: # if nan value
            worksheet.write(x+2,3,-1000,format_cell)
        worksheet.write(x+2,4,hydrophilicity_score[i],format_cell)
        worksheet.write(x+2,5,av_solvent_accessibility[i],format_cell)
        worksheet.write(x+2,6,av_score[i],format_cell)
        worksheet.write(x+2,7,round(pareto_dominant[i],2),format_cell)
        worksheet.write(x+2,8,max_repeat_length[i],format_cell)
        worksheet.write(x+2,9,i,format_cell)
        format_cell.set_font_color('black')
    workbook.close()
    
def plot_epitope_ranking(folder,project_name,av_score_isolates,hydrophilicity_score,av_score2,av_solvent_accessibility,mec):
    
    fig = plt.figure(num=None,figsize=(5,3),dpi=72)
    fig.add_axes([0.25,0.15,0.7,0.75])
    plt.title(project_name)
    XMIN = 5.2
    XMAX = 7.5        
    plt.xlim(XMIN,XMAX)
    plt.ylim(-1.2,1.5)
    plt.xlabel(r'Average pair-wise BL95 score within isolates')
    plt.ylabel(r'Hydrophilicity')

    plt.scatter(av_score_isolates,hydrophilicity_score,s=av_score2,c=av_solvent_accessibility,vmin=50,vmax=100.0,edgecolors=mec,cmap=plt.cm.get_cmap('gray'))
    
    #for kdx in range(nb_epitopes):
    #    plt.annotate(str(kdx),xy=(av_score_isolates[kdx],max_homology_identity[kdx]),xytext=(0,0),textcoords='offset points',color='k',fontsize=10) 
    plt.hlines(y=0,xmin=XMIN,xmax=XMAX,linewidth=0.5,linestyle='dotted',color='black') # 70% identity requirement
    cbar = plt.colorbar()
    cbar.set_label('Average solvent accessibility [%]')
    #plt.show()  
    fig.savefig(folder+'FIGURE-final-epitopes.png', dpi=300)

def get_indices_FIG(percent_conservation,indices_bs_variability,indices_H,indices_E,indices_DR,buried_score):
    
    L = len(percent_conservation)    
    indices_FIG = np.zeros([6,L],dtype=np.float)
    for i in range(L):
        if indices_H.count(str(i+1)) > 0:
            indices_FIG[0,i] = -4
        elif indices_E.count(str(i+1)) > 0:
            indices_FIG[0,i] = -6
        else:
            indices_FIG[0,i] = -5
            
        if indices_DR.count(str(i+1)) > 0:
            indices_FIG[2,i] = -3
        else:
            indices_FIG[2,i] = float('nan')
        ##########Modified by jibril################
        if i==len(buried_score):
            buried_score=np.append(buried_score,[np.nan])
        ############################################
        indices_FIG[1,i] = (buried_score[i]/3.0 - 2.5)
         
        if indices_bs_variability.count(str(i+1)) > 0: 
            indices_FIG[3,i] = 4.5 
        else:
            indices_FIG[3,i] = float('nan')
            
        indices_FIG[4,i] = float(percent_conservation[i])+2.5 # 3
            
    return indices_FIG
    
def perform_hierarchical_clustering_predictions(folder,methods_residues,dist_metric='jaccard',linkage_method='single'):
  
    # Hierarchical clustering
    DIST = ssd.squareform(ssd.pdist(methods_residues,dist_metric))
    Z0 = cl.linkage(ssd.squareform(DIST),linkage_method)
    plt.figure()
    d0 = cl.dendrogram(Z0)#,labels=ep.methods)
    #plt.show()
    print '..... Hierarchical clustering linkage matrix:\n'
    DIST = remove_diag(DIST,float('nan'))
    new_ordering = d0['leaves']
    new_labels = [ep.methods[i] for i in new_ordering]
    #with open(folder.split("/")[-3]+"_matrix.pickle", "w") as output_file:
    #    pickle.dump([DIST[:,new_ordering],new_labels], output_file)
    plot_matrix(folder,DIST[:,new_ordering],new_labels)
    
    return Z0
    
def compute_repeat_sequences(seq):

    size_repeat = 1
    aa = seq[0]
    size_repeat_v = []
    for i in range(1,len(seq)):
        if seq[i]==aa:
            size_repeat += 1
        else:
            aa = seq[i]
            size_repeat_v.append(size_repeat)
            size_repeat = 1
    size_repeat_v.append(size_repeat) # append as well the last value that was found 
            
    return size_repeat_v

#%% GLOBAL FUNCTIONS

def get_all_epitope_predictions(folder,data_folder,L,project_name):
    
    print('\n>> EPITOPE PREDICTION | RETRIEVE PREDICTIONS FROM ALL METHODS\n')
    
    epitopes = collections.defaultdict(dict) # compile results in a dictionary 
    
    for met in ep.methods:
        mymethod = {'method': eval(met+'_epitopes')}
        (ml,score_res) = mymethod['method'](data_folder,L)
        visualize_epitopes(ml,score_res,folder,project_name,met)
        epitopes[met] = {'method':met, 'ml':ml,'score_res':score_res}   
    #pickle.dump(epitopes, open(folder+'epitope_predictions.p','wb'))
    
    all_ml = []
    [ [all_ml.extend(epitopes[method]['ml'][j]) for j in range(len(epitopes[method]['ml']))] for method in ep.methods] 
    
    # majority score in [-6,6], then find regions with a positive score (i.e. residue is present as epitope in at least 3 methods)
    Nhalf = len(ep.methods)/2.0
    score_majority = []
    for kdx in range(L):
        indices_pos = [i for i, x in enumerate(all_ml) if x == kdx] # index of positive epitope fragments where the residue appears
        score_majority.append(len(indices_pos)) # positive scores of the epitope fragments where the residue appears
    score_majority = (np.array(score_majority)-Nhalf)/Nhalf
    score_majority_smooth = uf.smooth(score_majority.tolist(),5)
    
    input = (np.where(score_majority_smooth>0)[0]).tolist()
    ml_majority = uf.get_continuous_fragments(input)                     
    
    return (epitopes,ml_majority,score_majority)

def view_epitope_predictions(figname,folder,project_name,epitopes,ml_majority,score,SP,TM,PTM,percent_conservation,indices_bs_variability,indices_H,indices_E,indices_DR,buried_score,ylocs=[]):
    
    print('\n>> EPITOPE PREDICTION | VIEW PREDICTIONS FROM ALL METHODS\n')
    
    L = len(score)  
    N_methods = len(ep.methods)
    color_method = [list(plt.cm.jet(i)) for i in range(0,255,int(round(255/N_methods))) ]
    if ylocs == []:
        ylocs = range(125,125+(N_methods)*10,10)
        ylocs = [ylocs[i]/100.0 for i in range(len(ylocs))]

    #fig = plt.figure(num=None,figsize=(5,3),dpi=72)
    fig = plt.figure(num=None,figsize=(4,3),dpi=300)
    ax = fig.add_axes([0.25,0.15,0.7,0.75])
    plt.title(project_name)
    plt.xlim(0,L)
    plt.ylim(ep.YMIN,ep.YMAX)
        
    # for visualization            
    ax.add_patch(patches.Rectangle((0, 2.5),L,1,facecolor='lightgray',edgecolor='lightgray',alpha=0.5))
    ax.add_patch(patches.Rectangle((0, -2.5),L,1,facecolor='lightgray',edgecolor='lightgray',alpha=0.5))    
    plt.hlines(y=3,xmin=0,xmax=L,linewidth=0.5,linestyle='dotted',color='black')
    
    # signal peptide
    try:
        ax.add_patch(patches.Rectangle((SP[0], ep.YMIN),SP[-1]-SP[0], ep.YMAX-ep.YMIN,facecolor='lightgreen',edgecolor='lightgreen',alpha=0.5))
    except:
        pass
    
    # transmembrane domain
    try:
        ax.add_patch(patches.Rectangle((TM[0], ep.YMIN),TM[-1]-TM[0], ep.YMAX-ep.YMIN,facecolor='orange',edgecolor='orange',alpha=0.5))
    except:
        pass
            
    try:
        for kdx in PTM:
            plt.vlines(x=kdx,ymin=ep.YMIN,ymax=ep.YMAX,color='magenta')
    except:
        pass

    indices_FIG = get_indices_FIG(percent_conservation,indices_bs_variability,indices_H,indices_E,indices_DR,buried_score)

    plt.plot(indices_FIG[0,:],'k-')
    plt.plot(indices_FIG[1,:],'k.-',ms=1)
    plt.plot(indices_FIG[2,:],'k-')
    plt.plot(indices_FIG[3,:],'k-')
    plt.plot(indices_FIG[4,:],'k-')
    plt.plot(indices_FIG[5,:],'k-')
    plt.plot(score,'k-',linewidth=1,color='lightgray')
    plt.plot(uf.smooth(score,5).tolist(),'k-',linewidth=0.5)   
    plt.xticks(np.arange(0,L,50))
    ax.xaxis.grid() # vertical lines
    plt.xlabel(r'# amino acid')
    plt.yticks([-6,-5,-4,-3,-2,0,1.5,3,4.5])
    ax.set_yticklabels((r'$\beta$',r'coil',r'$\alpha$',r'DISORDERED',r'SOLV-ACC',r'RESIDUE-score',r'EPITOPE-pept',r'ISOLATES-CONS',r'HOMOL-VAR'))
    matplotlib.rc('font', size=5)
    if len(epitopes) > 0:    
        for i in range(N_methods):
            for kdx in range(len(epitopes[ep.methods[i]]['ml'])):
                plt.hlines(y=ylocs[i],xmin=epitopes[ep.methods[i]]['ml'][kdx][0],xmax=epitopes[ep.methods[i]]['ml'][kdx][-1],linewidth=.5,linestyle='solid',color=color_method[i])   
        for kdx in range(len(ml_majority)):
            ax.add_patch(patches.Rectangle((ml_majority[kdx][0],-1.1),ml_majority[kdx][-1]-ml_majority[kdx][0],3.2,facecolor='lightgray',edgecolor='lightgray',alpha=0.5))
    else:
        ax.add_patch(patches.Rectangle((0, -1),L,2,facecolor='lightgray',edgecolor='lightgray',alpha=0.5))
        for kdx in range(len(ml_majority)):
            loc = ylocs[kdx%len(ylocs)] 
            plt.hlines(y=loc,xmin=ml_majority[kdx][0],xmax=ml_majority[kdx][-1],linewidth=0.8,linestyle='solid',color='black')  
            #ax.add_patch(patches.Rectangle((ml_majority[kdx][0],-1.1),ml_majority[kdx][-1]-ml_majority[kdx][0],2.2,facecolor='lightgray',edgecolor='lightgray',alpha=0.5))
    plt.hlines(y=0,xmin=0,xmax=L,linewidth=0.5,linestyle='dotted',color='black')
    #plt.show()
    fig.savefig(folder+'FIGURE-'+figname+'.png', dpi=fig.dpi)
    
    output = 1
    return output

def get_consensus_epitopes(epitopes,L,score_majority,clustering_folder=[],folder=[]):
    
    print('\n>> EPITOPE PREDICTION | 2-BY-2 METHODS CORRELATION\n')   
    
    N_methods = len(ep.methods)
    
    # methods_residues contains the output of each prediction methods at every amino acid position (1,0)
    methods_residues = np.zeros((N_methods,L),dtype=int)
    for kdx in range(N_methods):        
        flat_epitopes = sum(epitopes[ep.methods[kdx]]['ml'],[])
        for i in range(L):
            if flat_epitopes.count(i)>0:
                methods_residues[kdx,i] = 1    
    
    # perform hierarchical clustering of all-methods predictions
#    if len(clustering_folder) > 0:
#        if hasattr(ep.DIST_METRIC, '__call__'): # i.e. ep.DIST_METRIC is a function
#            Z0 = np.load(clustering_folder+'methods_hierarchical_clustering_'+ep.DIST_METRIC.__name__+'_'+ep.LINKAGE_METHOD+'.npy')
#            print '..... Using '+ep.DIST_METRIC.__name__+' metric to get consensus'
#        else:
#            Z0 = np.load(clustering_folder+'methods_hierarchical_clustering_'+ep.DIST_METRIC+'_'+ep.LINKAGE_METHOD+'.npy')
#            print '..... Using '+ep.DIST_METRIC+' metric to get consensus'
#    else:
        Z0 = perform_hierarchical_clustering_predictions(folder,methods_residues,ep.DIST_METRIC,ep.LINKAGE_METHOD)

    # Use hierarchical tree to compute weighted predictions
    weighted_score_res = []  
    for pos in range(L):
        D = methods_residues[:,pos]
        weighted_score_res_tmp = []
        for i in range(np.shape(Z0)[0]):
            if (Z0[i,0] < N_methods) & (Z0[i,1] < N_methods):  # make a non-weighted average (0.5/0.5)
                weighted_score_res_tmp.append( (D[int(Z0[i,0])]+D[int(Z0[i,1])])/2.0 )
                
            elif (Z0[i,0] < N_methods):    
                Lx = Z0[i,2]
                mat_idx = int(Z0[i,1]%N_methods)
                x = float(Z0[mat_idx,2])
                if x>=Lx:
                    print '..... !!!!! PREDICT_EPITOPES.get_methods_correlation (a): x should always be smaller than L !!!!!'
                x = x/Lx*(2/3.0-0.5)+0.5  # map x from [0,L] to [0.5,2/3] 
                weighted_score_res_tmp.append(D[int(Z0[i,0])]*(1-x) + weighted_score_res_tmp[mat_idx]*x )
            else:
                #print '..... Linkage matrix: case with four arms .....'
                mat_idx0 = int(Z0[i,0]%N_methods)
                mat_idx1 = int(Z0[i,1]%N_methods)
                Lx = Z0[mat_idx1,2]
                x = float(Z0[mat_idx0,2])
                if x>=Lx:
                    print '..... !!!!! PREDICT_EPITOPES.get_methods_correlation (b): x should always be smaller than L !!!!!'
                x = x/Lx*(0.5-1/3.0)+1/3.0 # !! different re-scaling than from above is needed here !!
                weighted_score_res_tmp.append( weighted_score_res_tmp[mat_idx0]*x + weighted_score_res_tmp[mat_idx1]*(1-x) )
        weighted_score_res.append(weighted_score_res_tmp[-1])
    
    weighted_score_res = np.array(weighted_score_res)
    weighted_score_res = weighted_score_res*2.0 - 1.0   # map to [-1,1]
    
#    plt.figure()
#    plt.plot(weighted_score_res)
#    plt.plot(score_majority,'r-')
#    plt.xlim(0,L)
#    plt.xlabel('#amino acid')
#    plt.ylabel('weigthed (b) & majority (r) score')
#    plt.hlines(y=0,xmin=0,xmax=L)
    
    return weighted_score_res

def crop_epitopes(score_res):

    print('\n>> EPITOPE PREDICTION | CROP EPITOPES\n')   

    lmax = ep.max_epitope_length
    lmin = ep.min_epitope_length

    score = uf.smooth(score_res,3).tolist()
    L = len(score)        
        
    ml_vec = []
    peak = np.max(score) # arbitrary positive value to start
    frag = []
    frag.append(np.argmax(score))
    score_frag = peak

    # Extend the peak
    while peak > 0: 
        while ((len(frag)<lmax) & (score_frag>0)):    
    
            if ((frag[0]>0) & (frag[-1]<L-1)): # check that we haven't reached the first or last amino acid of the protein, in which case there is only one left direction to go for expanding
                current = [frag[0]-1,frag[-1]+1]
            elif ((frag[0]==0) & (frag[-1]<L-1)):
                current = [frag[-1]+1,frag[-1]+1]
            else:
                current = [frag[0]-1,frag[0]-1]
                
            i_tmp = np.argmax([score[current[0]],score[current[1]]]) # = 0 or 1 (i.e. expand to left (0) or right (1))
            frag.append(current[i_tmp])
            frag.sort()  
            score_frag = score[current[i_tmp]]
            
        ml_vec.append(frag)
        
        score2 = np.array(score)
        score2[list(set(sum(ml_vec,[])))] = -2 # set the fragment to -2 to find the new peak
        
        peak = np.max(score2)
        idx_peak = np.argmax(score2)
        
        frag = []
        frag.append(idx_peak)
        score_frag = peak
        
    # Keep only those fragments whose length is between min_length-max_length
    ml_vec_long = [ml_vec[i] for i in range(len(ml_vec)) if len(ml_vec[i])>=lmin]   
    
    # Control: check that the maximal length is 20aa
    ml_vec_L = [len(ml_vec_long[i]) for i in range(len(ml_vec_long))]
    if (max(ml_vec_L) > lmax) | (min(ml_vec_L) < lmin):
        print '..... !!!!! PREDICT_EPITOPES.crop_epitopes: something weird !!!!!'
     
    # sort epitopes by score, for better visualization 
    score = np.array(score)        
    av_score = [ np.mean(score[ml_vec_long[i]]) for i in range(len(ml_vec_long))]
    idx_sort = np.argsort(av_score)[::-1] # sort in ascending order and then reverse the indices in the list (by using [::-1]); sorting will be useful for visualization
    ml_vec_long = np.array(ml_vec_long)[idx_sort].tolist() # sort ml_vec_long accoring to score (top scoring on top)          
    av_score = np.sort(av_score)[::-1].tolist() 
        
    return ml_vec_long
    
def post_process_cysteins(ml,refseq,Hphob_scale):
    
    print('\n>> EPITOPE PREDICTION | CYSTEINS POST-PROCESSING\n')   
    
    PERC = ep.PERC # percentage that defines the "border" of the epitope, where we would like to identify cysteins to postprocess the epitope
    
    nb_epitopes = len(ml)    
    
    cystein_status = []                                                     # 0 if no cystein, -1 if cystein in the middle or post-processing would shorten too much the epitope, 1 if cystein could be processed
    ml_processed = []        
    for kdx in range(nb_epitopes):
        epitope_sequence = refseq[ml[kdx][0]:ml[kdx][-1]+1]
        nb_cysteins = epitope_sequence.count('C')
        if nb_cysteins > 0: # check where they are
            pos_cysteins = [i for i in range(len(epitope_sequence)) if epitope_sequence[i]=='C']
            pos_middle_cysteins = list(set(pos_cysteins) & set(range(int(PERC*len(epitope_sequence)),len(epitope_sequence)-int(PERC*len(epitope_sequence))))) # check if one of the cysteins is in the middle (not good)
            pos_Nter_cysteins = list(set(pos_cysteins) & set(range(0,int(PERC*len(epitope_sequence))))) 
            pos_Cter_cysteins = list(set(pos_cysteins) & set(range(len(epitope_sequence)-int(PERC*len(epitope_sequence)),len(epitope_sequence)))) 
            if len(pos_middle_cysteins) > 0:                                # i.e. cysteins in the middle of the epitope
                cystein_status.append(-1)
                ml_processed.append(ml[kdx])
            else:                                                           # try post-processing (check length is still above 15aa after post-processing)
            
                if len(pos_Nter_cysteins)==nb_cysteins:                     # all cysteins are close to the N-terminus
                    new_start = ml[kdx][max(pos_cysteins)]                  # new start position for that epitope
                    new_end = ml[kdx][-1]                                    
                    for idxaa in range(max(pos_cysteins)):                  # try extending epitope on the C-ter side             
                        if float(Hphob_scale[refseq[ml[kdx][-1]+1+idxaa]]) <= 0:  # i.e. neutral or hydrophilic aa
                            new_end = ml[kdx][-1]+1+idxaa
                        else:                                               # break to stop adding amino acids
                            break
                elif len(pos_Cter_cysteins)==nb_cysteins:                   # all cysteins are close to the C-terminus
                    new_start = ml[kdx][0]
                    new_end = ml[kdx][min(pos_cysteins)]
                    for idxaa in range(len(epitope_sequence)-min(pos_cysteins)-1):                          
                        if float(Hphob_scale[refseq[ml[kdx][0]-1-idxaa]]) <= 0:
                            new_start = ml[kdx][0]-1-idxaa
                        else:
                            break
                else:                                                       # there are cysteins on both termini
                    new_start = max(pos_Nter_cysteins)
                    new_end = min(pos_Cter_cysteins)
                    
                if new_end-new_start+1 >= ep.min_epitope_length:             # check that processed epitopes are long enough
                    cystein_status.append(1)
                    old_score = [Hphob_scale[refseq[ml[kdx][i]]] for i in range(len(ml[kdx]))]
                    #print '..... '+str(ml[kdx][0])+' --> '+str(new_start)
                    #print '..... '+str(ml[kdx][-1])+' --> '+str(new_end)
                    ml_processed.append(np.array(range(new_start,new_end+1)))
                    new_score = [Hphob_scale[refseq[ml_processed[kdx][i]]] for i in range(len(ml_processed[kdx]))]
                    #print '..... '+str(np.mean(old_score))+' --> '+str(np.mean(new_score))+'\n'
                else:
                    cystein_status.append(-1)
                    ml_processed.append(ml[kdx])
        else:                                                               # i.e. there are no cysteins in the epitope sequence
            cystein_status.append(0)
            ml_processed.append(ml[kdx])
    
    return ml_processed

def terminal_epitopes(folder,ml_processed,seq_ref,SP,Hphob_scale):

    print('\n>> EPITOPE PREDICTION | TERMINAL EPITOPES\n')  
    
    L = len(seq_ref)
    L_TR = ep.TERMINAL_REGION
    L_TR_MIN = ep.MIN_TERMINAL_REGION
    
    if SP==[]:
        SP = [-1]

    N_ter = seq_ref[SP[-1]+1:SP[-1]+L_TR+1]
    C_ter = seq_ref[-L_TR:]
    
    N_ter_kd = [float(Hphob_scale[N_ter[i]]) for i in range(len(N_ter))]
    C_ter_kd = [float(Hphob_scale[C_ter[i]]) for i in range(len(C_ter))]

    fig = plt.figure(num=None,figsize=(5,3),dpi=72)
    fig.add_axes([0.25,0.15,0.7,0.75])
    plt.xlim([SP[-1]+1,SP[-1]+21])
    plt.title('N-ter')
    plt.xlabel('# amino acid')
    plt.ylabel('Hydrophobicity score')
    plt.plot(range(SP[-1]+1,SP[-1]+L_TR+1),N_ter_kd,'b-')
    plt.hlines(y=0,xmin=1,xmax=L_TR,linewidth=0.5,linestyle='dotted',color='blue')
    #plt.show()
    fig.savefig(folder+'N-terminal-Hphob-score.pdf', dpi=fig.dpi)
    
    fig = plt.figure(num=None,figsize=(5,3),dpi=72)
    fig.add_axes([0.25,0.15,0.7,0.75])
    plt.xlim([L-L_TR,L])
    plt.title('C-ter')
    plt.xlabel('# amino acid')
    plt.ylabel('Hydrophobicity score')
    plt.plot(range(L-L_TR,L),C_ter_kd,'b-')
    plt.hlines(y=0,xmin=L-L_TR,xmax=L,linewidth=0.5,linestyle='dotted',color='blue')
    #plt.show()
    fig.savefig(folder+'C-terminal-Hphob-score.pdf', dpi=fig.dpi)
    
    H_score_N = []
    H_score_C = []
    ml_w_N = []
    ml_w_C = []
    windows = range(L_TR_MIN,L_TR+1)
    for kdx in windows:
        ml_tmp_N = range(SP[-1]+1,SP[-1]+1+kdx)
        ml_w_N.append(ml_tmp_N)
        
        ml_tmp_C = range(windows[-1]-kdx+(L-L_TR),windows[-1]+(L-L_TR))
        ml_w_C.append(ml_tmp_C)
        
        tmp_score = [float(Hphob_scale[N_ter[i-ml_tmp_N[0]]]) for i in ml_tmp_N]
        H_score_N.append(np.mean(tmp_score))
        
        if np.mean(tmp_score) < 0:
            ml_processed.append(ml_tmp_N)
        
        tmp_score = [float(Hphob_scale[C_ter[i-(L-L_TR)]]) for i in ml_tmp_C]
        H_score_C.append(np.mean(tmp_score))
        
        if np.mean(tmp_score) < 0:
            ml_processed.append(ml_tmp_C)
    
    fig = plt.figure(num=None,figsize=(5,3),dpi=72)
    fig.add_axes([0.25,0.15,0.7,0.75])
    plt.xlim([SP[-1]+1,SP[-1]+20])
    plt.title('N-ter')
    plt.xlabel('# amino acid')
    plt.ylabel('Hydrophobicity score')
    for r in range(len(ml_w_N)):
        plt.hlines(y=H_score_N[r],xmin=ml_w_N[r][0],xmax=ml_w_N[r][-1],color='blue')
    #plt.show()
    fig.savefig(folder+'N-terminal-epitopes.pdf', dpi=fig.dpi)
        
    fig = plt.figure(num=None,figsize=(5,3),dpi=72)
    fig.add_axes([0.25,0.15,0.7,0.75])
    plt.xlim([L-L_TR,L-1])
    plt.title('C-ter')
    plt.xlabel('# amino acid')
    plt.ylabel('Hydrophobicity score')
    for r in range(len(ml_w_C)):
        plt.hlines(y=H_score_C[r],xmin=ml_w_C[r][0],xmax=ml_w_C[r][-1],color='blue')
    #plt.show()
    fig.savefig(folder+'C-terminal-epitopes.pdf', dpi=fig.dpi)
    
    return ml_processed
    
def discard_epitopes_with_repeats(ml,refseq):
    
    nb_epitopes = len(ml)    
    
    # calculate repeat sequences
    max_repeat_length = [max(compute_repeat_sequences(refseq[ml[i][0]:ml[i][-1]+1])) for i in range(nb_epitopes)]
    
    ml_norepeats = [ml[i] for i in range(nb_epitopes) if max_repeat_length[i]<=ep.MAX_REPEAT_LENGTH]
    max_repeat_length = [max_repeat_length[i] for i in range(nb_epitopes) if max_repeat_length[i]<=ep.MAX_REPEAT_LENGTH] # remove same elements in max_repeat_length
    
    return (ml_norepeats,max_repeat_length)
    
################Modified by jibril
#def write_epitope_ws_fasta(folder,list_seq,uniprotID,ml_processed,genotype,extra_info):
#    
#    for i in range(len(ml_processed)):
#        fo = open(folder+'epitope-'+str(i)+'.fasta','w')
#        for j in range(len(list_seq)):
#            fo.write('>'+uniprotID[j]+'|'+genotype[j]+extra_info[j]+'\n')
#            seq = ''.join((list_seq[j,:]).tolist())
#            fo.write(seq[ml_processed[i][0]:ml_processed[i][-1]+1]+'\n')
#        fo.close()

#    output = 1
#    return output 

def write_epitope_ws_fasta(folder,list_seq,uniprotID,ml_processed):
    for i in range(len(ml_processed)):
        fo = open(folder+'epitope-'+str(i)+'.fasta','w')
        for j in range(len(uniprotID)):
            fo.write('>'+uniprotID[j]+'|'+'\n')
            seq = ''.join((list_seq[j,:]).tolist())
            fo.write(seq[ml_processed[i][0]:ml_processed[i][-1]+1]+'\n')
        fo.close()

    output = 1
    return output 


def draw_epitope_ws_weblogos(folder,nb_epitopes):
    
    print('\n>> EPITOPE PREDICTION | DRAWING WEBLOGOS\n')   

    for i in range(nb_epitopes):       
        fo = open(folder+'epitope-'+str(i)+'.fasta','r')
        #os.system(ep.seqlogo_cmd+' -k 0 -c -w 15 -h 4 -f '+folder+'epitope-'+str(i)+'.fasta > '+folder+'epitope-'+str(i)+'.ps')  
        os.system(ep.seqlogo_cmd+' -A "protein" -W 15 -f '+folder+'epitope-'+str(i)+'.fasta > '+folder+'epitope-'+str(i)+'.ps')  
        os.system('ps2pdf -dEPSCrop '+folder+'epitope-'+str(i)+'.ps '+folder+'epitope-'+str(i)+'.pdf')
        fo.close()
    
    output = 1
    return output
    
def blast_all_epitopes(folder,data_folder,ml,refseq):
    
    print('\n>> EPITOPE PREDICTION | BLAST ALL EPITOPES (be patient...)\n')   
    
    fo = open(folder+'REF-SEQ.fasta','w')
    fo.write('>refseq\n'+refseq)
    fo.close()    
    
    nb_epitopes = len(ml)
#    max_homology_identity = []
    for kdx in range(nb_epitopes): # perform a BLAST for each epitope fragment
#        fname = folder+'epitope'+str(kdx)+'_BLAST.xml'
        fname_TXT = folder+'epitope'+str(kdx)+'_BLAST.txt'
        print "-----------BLAST ",fname_TXT
#        if not os.path.isfile(fname): # do not run BLAST if already ran once
#            print '..... Performing BLAST for epitope # '+str(kdx)+'/'+str(nb_epitopes-1)+'... (xml)'
#            blast_cline = NcbiblastpCommandline(cmd=ep.blast_loc,query=folder+'REF-SEQ.fasta',remote=True,out=fname,import_search_strategy=data_folder+'blast_search_strategy.asn',query_loc=str(ml[kdx][0]+1)+'-'+str(ml[kdx][-1]+1),outfmt=5)
#            blast_cline()
        if not os.path.isfile(fname_TXT): # do not run BLAST if already ran once
            print '..... Performing BLAST for epitope # '+str(kdx)+'/'+str(nb_epitopes-1)+'... (txt)'
            blast_cline = NcbiblastpCommandline(cmd=ep.blast_loc,query=folder+'REF-SEQ.fasta',remote=True,out=fname_TXT,import_search_strategy=data_folder+'blast_search_strategy.asn',query_loc=str(ml[kdx][0]+1)+'-'+str(ml[kdx][-1]+1),outfmt=0)
            blast_cline()
#        blast_qresult = SearchIO.read(fname, 'blast-xml')
#        identity = [100.0*blast_qresult[i][0].ident_num for i in range(len(blast_qresult))]
        #print max(np.array(identity)/float(ml[kdx][-1]-ml[kdx][0]+1))
#        if len(identity)>0:
#            max_homology_identity.append(max(identity)/float(ml[kdx][-1]-ml[kdx][0]+1))
#        else:
#            max_homology_identity.append(0)
            
#    return max_homology_identity
    
def rank_epitopes(folder,refseq,project_name,ml,score,av_score_conservation,buried_score,hydrophilicity_score_pos,max_repeat_length,indices_SIGNAL,indices_TM,indices_PTM):
    
    print('\n>> EPITOPE PREDICTION | RANK EPITOPES\n')   
    
    nb_epitopes = len(ml) 
    
    #score = (score+1)/2.0 # rescale [-1;1] to [0;1]
    av_score = [np.mean(score[ml[i]]) for i in range(nb_epitopes)]
    av_score2 = [(av_score[i]+1)/2.0 for i in range(nb_epitopes)] 
    av_score2 = [int(round(av_score2[i]*40 + 4)*5) for i in range(nb_epitopes)] # rescale [0,1] to marker size of [4, 24] i.e. increments of 0.1
    
    av_score_isolates = [np.nanmean(np.array(av_score_conservation)[ml[i]]) for i in range(nb_epitopes)]
    
    buried_score[np.where(buried_score==1)] = 5
    buried_score[np.where(buried_score==2)] = 25
    buried_score[np.where(np.isnan(buried_score))] = 100
    av_solvent_accessibility = [int(round(np.mean(buried_score[ml[i]]))) for i in range(nb_epitopes)]

    hydrophilicity_score = [np.mean(hydrophilicity_score_pos[ml[i]]) for i in range(nb_epitopes)]

    # Check if residue overlaps with signal peptide, transmembrane domain or harbours post-translational modifications
    mec = check_overlap(ml,indices_SIGNAL,indices_TM,indices_PTM)
    # Compute Pareto dominance
    data = zip(av_score_isolates,hydrophilicity_score,av_score,av_solvent_accessibility) # goal: maximize each of these features
    (idx_pareto,pareto_dominant) = compute_pareto_dominance(data)
    
    # Print epitope ranking to excel file
    write_excel_file(folder,project_name,ml,idx_pareto,mec,refseq,av_score_isolates,hydrophilicity_score,av_score,av_solvent_accessibility,pareto_dominant,max_repeat_length)
    
    # Sort by largest score before plotting (so big markers do not hide smaller ones below them)         
    idx_sort = np.argsort(av_score)[::-1]
    av_score = np.array(av_score)[idx_sort].tolist()
    av_score_isolates = np.array(av_score_isolates)[idx_sort].tolist()
    hydrophilicity_score = np.array(hydrophilicity_score)[idx_sort].tolist()
    av_solvent_accessibility = np.array(av_solvent_accessibility)[idx_sort].tolist()
    mec = np.array(mec)[idx_sort].tolist()
    max_repeat_length = np.array(max_repeat_length)[idx_sort].tolist()
    ml = np.array(ml)[idx_sort].tolist() # re-sort ml as well    
    
    # Plot result and save figure
    plot_epitope_ranking(folder,project_name,av_score_isolates,hydrophilicity_score,av_score2,av_solvent_accessibility,mec)
    
    return (data,idx_pareto,pareto_dominant)
    
def get_all_sliding_windows(refseq,SP):
    
    sliding_windows = []    

    # add sliding at the terminal regions, which are allowed to be shorter    
    L = len(refseq)
    L_TR = ep.TERMINAL_REGION
    L_TR_MIN = ep.MIN_TERMINAL_REGION
    
    if SP==[]:
        SP = [-1]
    
    windows = range(L_TR_MIN,L_TR+1)
    for kdx in windows:
        ml_tmp_N = range(SP[-1]+1,SP[-1]+1+kdx)
        sliding_windows.append(ml_tmp_N)
        
        ml_tmp_C = range(windows[-1]-kdx+(L-L_TR),windows[-1]+(L-L_TR))
        sliding_windows.append(ml_tmp_C)
        
    # add sliding windows of max_epitope_length for the rest of the protein
    for kdx in range(SP[-1]+2,len(refseq)-ep.max_epitope_length):
        sliding_windows.append(range(kdx,kdx+ep.max_epitope_length))
    
    return sliding_windows
    
def keep_epitopes_with_conserved_hydrophilic_stretches(ml,percent_conservation,thres,hphil_score):
    
    ml_long_stretches = []
    pos_cons_stretches = []
    hphil_cons_stretches= []
    L_stretch = []    
    
    percc = np.array(percent_conservation)
    hphil_score = np.array(hphil_score)
    for w in ml:
        w_cons = (percc[w]>=thres).astype(int)
        w_cons = [i for i,x in enumerate(w_cons) if x==1]
        ml_all = uf.get_continuous_fragments(w_cons)
        if len(ml_all)>0:
            ml_long_L = max([len(mla) for mla in ml_all]) # longest stretch of consecutively conserved aa
            ml_long = [mla for mla in ml_all if len(mla)==ml_long_L][0] # longest stretch of consecutively conserved aa
            av_hphil_score = np.mean(hphil_score[np.array(w)[ml_long]]) # average hydrophilicity of longest stretch of consecutively conserved aa      
            if (ml_long_L >= ep.MIN_STRETCH_LENGTH and av_hphil_score>=0):
                ml_long_stretches.append(w)
                pos_cons_stretches.append(np.array(w)[ml_long])
                hphil_cons_stretches.append(av_hphil_score)
                L_stretch.append(ml_long_L)
    return (ml_long_stretches,pos_cons_stretches,hphil_cons_stretches,L_stretch)
    
def discard_epitopes_with_PTM(ml,indices_PTM):
    
    ml_noPTM = []
    
    for mli in ml:
        ptm = 0
        for i in indices_PTM:
            if mli.count(i)>0:
                ptm = 1
                break
        if ptm==0:
            ml_noPTM.append(mli)
    
    return ml_noPTM
        
def discard_hydrophobic_epitopes(ml,hphil_score):
    
    ml_hphile = []   
    
    for mli in ml:
        av_hphil_score = np.mean( hphil_score[mli] )  
        if av_hphil_score>0:
            ml_hphile.append(mli)
    
    return ml_hphile
