# -*- coding: utf-8 -*-
"""
@author: aitana neves | heig-vd
"""

global retrieve_homologous_sequences
global align_all_sequences
global compute_position_specific_variability

#%% IMPORTS

import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import epitope_prediction as ep
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio.SubsMat import MatrixInfo as sm
from epitope_prediction import useful_functions as uf
import numpy as np
import time

#%% GLOBAL FUNCTIONS

def retrieve_homologous_sequences(folder,refseq,myorganism):

    print('\n>> BETWEEN-SPECIES | RETRIEVING BLAST HOMOLOGOUS SEQUENCES\n')
    myorganism=myorganism.replace("+"," ")
    print myorganism
    fname = folder+'homologous_sequences.xml'
    
   # if not os.path.isfile(fname): # do not run BLAST if already run once  
    aln=[[]]
    stop=True
    while stop:
        # if not os.path.isfile(fname): # do not run BLAST if already run once
        try :
            aln=[[]]
            if len(myorganism.split(" "))<1:
                result_blast = NCBIWWW.qblast('blastp', 'refseq_protein', refseq.strip('\n'),matrix_name='BLOSUM80',expect=0.001)
            else:
                result_blast = NCBIWWW.qblast('blastp', 'refseq_protein', refseq.strip('\n'),matrix_name='BLOSUM80',expect=0.001,entrez_query='NOT '+myorganism+'[organism]')
            print '..... Blast done.'        
            fo = open(fname,'w')        
            fo.write(result_blast.read()) 
            fo.close()
            result_blast.close()
            #else:
             #   print '..... Blast had already been performed. Retrieving existing results.'
            result_blast = open(fname)
            blast_records = NCBIXML.parse(result_blast)
        
            aln = [blast_record.alignments for blast_record in blast_records]
            if len(aln[0])>2 or len(myorganism.split(" "))<1:
                stop=False
            myorganism=" ".join(myorganism.split(" ")[:-1])
        except :
            time.sleep(3)
            stop = True
            fo = open(fname,'w')
            fo.close()
            result_blast.close()
    
    gis = [str(aln[0][i].title).split('|')[1] for i in range(len(aln[0]))]
    # get whole genomic sequence (not just where hit) for each hit (given by its accession number)
    Entrez.email = 'jibril.mammeri@heig-vd.ch'
    request = Entrez.epost("protein",id=",".join(map(str,gis)))
    result = Entrez.read(request)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    handle = Entrez.efetch(db="protein",retmode="xml", webenv=webEnv, query_key=queryKey)
    seq_blast = []
    for r in Entrez.parse(handle):
        # Grab the GI 
        try:
            gi=int([x for x in r['GBSeq_other-seqids'] if "gi" in x][0].split("|")[1])
        except ValueError:
            gi=None
            time.sleep(1)
        seq_blast.append('>gi|'+str(gi)+' | '+r['GBSeq_definition']+'\n'+r['GBSeq_sequence'].upper())
    seq_blast[:0] = ['>refseq\n'+refseq.strip('\n')]  # add query "ref" sequence
    tmp = "\n".join(seq_blast)
    fo = open(fname[0:-3]+'fasta','w') # write a fasta file
    fo.write(tmp)   
    fo.close()
    
    return fname
 
def align_all_sequences(fname): 
    stderr=None
    if not os.path.isfile(fname[0:-3]+'aln'):
        print '\n-- Performing clustalw alignment 2 --\n'
        clustalw_cline = ClustalwCommandline(cmd=ep.clustalw_loc, infile=fname[0:-3]+'fasta', outorder='input')
        stdout, stderr = clustalw_cline()   
    clustalw_align = AlignIO.read(fname[0:-3]+'aln', 'clustal')
    print('..... Alignment:\n')    
    #print(clustalw_align) 
    
    N_aln = np.shape(clustalw_align)[1]        
    list_homo_seq = np.zeros((len(clustalw_align),N_aln),dtype=np.string0)    
    list_homo_id = []
    for record in range(len(clustalw_align)): # loop over each nr sequence
        list_homo_id.append(clustalw_align[record].description)
        list_homo_seq[record,:] = (list(clustalw_align[record].seq))         
         
    indices_homo = [i for i, x in enumerate((list_homo_seq[0,:]).tolist()) if x != '-'] # Find not '-' in the reference sequence  

    return (list_homo_seq,indices_homo,stderr)          
  
def compute_position_specific_variability(folder,list_homo_seq,indices_homo):

    BL = sm.blosum95
    
    residue_conserved = [] # 1 if position consists of >95% likely mutations, 0 otherwise
    percent_conservation = [] # percent of conserved residues at each position
    score_conservation = []  # average pair-wise blosum score at each position (for the pairs all-vs-ref)  
  
    residue_conserved = []
    neg_scores= []
    for idxL in indices_homo: # loop over amino acids where there is no '-' in the reference sequence (i.e. in VHSV or IHNV)
   
        # determine score of each position using BLOSUM95 matrix    
        SEQ_homo = str((list_homo_seq[:,idxL]).tolist()).replace("'","").replace(",","").replace(' ','').strip('[]')     # positions were reference sequences has no '-' (useful for visualization in Protter)      
        score_tmp = []
        for jdxL2 in range(len(SEQ_homo)): # pair-wise computation of BLOSUM score against the reference genome
            try:
                score_tmp.append(BL[SEQ_homo[1],SEQ_homo[jdxL2]]) # use ref sequence in VHSV or IHNV and compare against it (i.e. not for all possible pairs)
            except:
                score_tmp.append(-8) # penalty for a gap in another species
        score_conservation = np.mean(score_tmp)        
        neg_scores = [i for i, x in enumerate(score_tmp) if x < 0] 
        percent_conservation.append(float(len(neg_scores))/np.shape(list_homo_seq)[0])             
        if score_conservation < 0: # if average score is negative (better than using 5%, because number of sequences is between 4-12, so 5% is less than one sequence)
            residue_conserved.append(1)
        else:
            residue_conserved.append(0)     
    
    indices_bs_variability = uf.indices2string([i for i, x in enumerate(residue_conserved) if x > 0]) # == indices_VAR_REFSEQ       
    indices_bs_variability = indices_bs_variability.split(',')    

    return (residue_conserved,percent_conservation,indices_bs_variability)
