# -*- coding: utf-8 -*-
"""
@author: aitana neves | heig-vd
"""

from epitope_prediction import useful_functions as uf

#%% EPITOPE_PREDICTION package
#----------------------------------------------------------------------------

# subfolders to be created
sf = ['within-species_protein_sequences','between-species_protein_sequences','all-methods-epitope-predictions','terminal-regions-epitopes','weblogos_epitopes','blast_epitopes','final_epitope_predictions_ranking']
sf2 = ['within-species','between-species','epitope-predictions-from-online-methods','weblogos-epitopes','final-predictions']

# commands
muscle_loc = r'/usr/bin/muscle'
clustalw_loc = r'/usr/bin/clustalw'
blast_loc = r'blastp'
seqlogo_cmd = r'weblogo'

# required threshold of conservation to determine conserved stretches in sliding windows (main2)
#THRES_CONSERVATION = 1.0

# epitope prediction methods to be used in the analysis
methods = ['BepiPred','ABCPred','AAP','BCPred','SVMTriP']
#methods = ['BepiPred','ABCPred','AAP','BCPred','Epitopia','COBEPro','TEPRF','CBTope','SVMTriP']
#methods = ['BepiPred','ABCPred','AAP','BCPred','Epitopia','COBEPro','TEPRF','CBTope']#,'SVMTriP']
#methods = ['BepiPred','ABCPred','AAP','BCPred','COBEPro','TEPRF','CBTope','SVMTriP']

# min, max epitope length
#min_epitope_length = 15
#max_epitope_length = 20

MAX_REPEAT_LENGTH = 2 # max length of repeats in epitopes

#MIN_STRETCH_LENGTH = 12 # minimum number of consecutive conserved regions to keep a window as a potential epitope

# hydrophobicity scales
cowan_whittaker_Hphob_scale = {'A':0.35, \
                            'R':-1.5, \
                            'N':-0.99, \
                            'D':-2.150, \
                            'C':0.760, \
                            'E':-0.93, \
                            'Q':-1.95, \
                            'G':0, \
                            'H':-0.65, \
                            'I':1.83, \
                            'L':1.8, \
                            'K':-1.54, \
                            'M':1.1, \
                            'F':1.69, \
                            'P':0.84, \
                            'S':-0.63, \
                            'T':-0.27, \
                            'W':1.35, \
                            'Y':0.39, \
                            'V':1.32}

# the scale was re-centered (i.e. minus the average)
roseman_Hphob_scale = {'A':0.390, \
                            'R':-3.950, \
                            'N':-1.910, \
                            'D':-3.810, \
                            'C':0.250, \
                            'E':-1.300, \
                            'Q':-2.910, \
                            'G':0.000, \
                            'H':-0.640, \
                            'I':1.820, \
                            'L':1.820, \
                            'K':-2.770, \
                            'M':0.960, \
                            'F':2.270, \
                            'P':0.990, \
                            'S':-1.240, \
                            'T':-1.000, \
                            'W':2.130, \
                            'Y':1.470, \
                            'V':1.300}
                            
PERC = 0.2 # percentage that defines the "border" of the epitope, where we would like to identify cysteins to postprocess the epitope. predict_epitoeps.post_process_cysteins() 

CBTope_THRESHOLD = 4 # predict_epitopes.CBTope_epitopes()

#TERMINAL_REGION = 20 # length of the terminal regions (predict_epitopes.terminal_epitopes())
#MIN_TERMINAL_REGION = 12

YMIN = -6.5 # for epitope visualization plots. predict_epitopes.view_all_epitope_predictions()
YMAX = 5.5

DIST_METRIC = uf.jaccard_shoulder #'jaccard' # to compare epitope prediction methods. E.g. 'jaccard' or 'jaccard_shoulder'
LINKAGE_METHOD = 'complete' # single, complete