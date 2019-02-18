# -*- coding: utf-8 -*-
"""
@author: aitana neves | heig-vd
"""

global retrieve_all_sequences
global keep_only_subgenotypes_from_fishpathogens
global align_all_sequences
global compute_position_specific_variability
global get_ws_output
global get_conserved_positions

#%% IMPORTS

#import code # code.interact(local=locals())
import re
import numpy as np   
import urllib
from xml.dom import minidom

import epitope_prediction as ep
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

from Bio.SubsMat import MatrixInfo as sm
from epitope_prediction import useful_functions as uf
import matplotlib.pyplot as plt

from Bio import Entrez
#%% LOCAL FUNCTIONS

def process_url(name):
    
    name = name.lower()
    name = name.strip(' ')
    if name.count('(')>0:
        name = name[0:name.find('(')-1]
    name = name.replace(' ','+')
    return name

def retrieve_gene_name_N_organism(uniprot_ref):
    
    print('\n>> WITHIN-SPECIES | RETRIEVING UNIPROT PROTEIN SEQUENCES FOR '+str(uniprot_ref)+'\n')    
    
    url = 'http://www.uniprot.org/uniprot/'+uniprot_ref+'.xml'
    dom = minidom.parse(urllib.urlopen(url))
    c=0
    for node in dom.getElementsByTagName('gene'):
       c+=1
       for child in node.getElementsByTagName('name'):
            gene = str(child.firstChild.nodeValue)
    for node in dom.getElementsByTagName('protein'):
       for child in node.getElementsByTagName('recommendedName'):
           for grandchild in child.getElementsByTagName('fullName'):
               protein = str(grandchild.firstChild.nodeValue)       
               
    proteins = protein.strip(' ').split(' ')
    protein = 'name:'+proteins[0]
    for i in range(1,len(proteins)):
        protein = protein+r'+OR+name:'+proteins[i]
    if c>0:
        myname = protein+'+OR+gene:'+gene
    else:
        myname = protein
    #myname = 'gene:'+gene
    print('..... Search options : '+myname.replace('+',' ')+'\n')  
      
    myorganism = []        
    for node in dom.getElementsByTagName('organism'):
       for child in node.getElementsByTagName('name'):
            myorganism.append(child.firstChild.nodeValue)
    TaxID=[]
    for node in dom.getElementsByTagName('organism'):
        for child in node.getElementsByTagName('dbReference'):
            TaxID.append(child.getAttribute("id"))
       #for child in node.getElementsByTagName('dbReference'):
       #     myorganism.append(child.attributes['id'].value)
    myorganism = process_url(myorganism[0]) 
    taxid=str(TaxID[0])
    print('..... Search options : '+myorganism.replace('+',' ')+'\n')  
     
    return (myname,myorganism,taxid)

     
def write_file(myfile,data,dim=1):

    fo = open(myfile,'w') 
    extension = myfile.split('.')[-1]

    if dim==1:
        if extension == 'txt':        
            fo.write(str(data).strip('[]').replace(', ','\n').replace("'",''))
        elif extension == 'fasta': 
            fo.write('>'+str(data).strip('[]').replace('\\n','\n').replace(', ','>').replace("'",''))
            
    elif dim==2:
        [fo.write(str(data[i]).replace(', ',r'\t').replace("'",'').strip('[]')+r'\n') for i in range(len(data))]
    else:
        print '----- !!!!! within_species.write_file() : dim should be equal to 1 or 2 !!!!! -----'
        
    fo.close() 
          
def rewrite_file(myfile,is_in,genotype=[]):
    
    print '..... Rewriting files...\n'

    fo = open(myfile,'r')      
    extension = myfile.split('.')[-1]

    if extension == 'txt':    
        data = np.array(fo.read().splitlines())[is_in].tolist()
        fo.close()
        fow = open(myfile,'w') 
        fow.write(str(data).strip('[]').replace(', ','\n').replace("'",''))
    elif extension == 'fasta':
        data = (np.array(fo.read().split('>')))[1:] # because of split('>'), the first element is empty
        data = [str(data[is_in[i]][0:10])+str(genotype[i])+' | '+str(data[is_in[i]][10:]) for i in range(len(is_in))]
        #data = (data[is_in]).tolist()     
        fo.close()
        fow = open(myfile,'w')
        fow.write('>'+str(data).strip('[]').replace('\\n','\n').replace(', ','>').replace("'",''))
    fow.close()
    
    return data

#%% GLOBAL FUNCTIONS

def retrieve_all_sequences_uniprot(folder,uniprot_ref,myname,myorganism,refseq):
    (myname_tmp,organism_tmp,taxid) = retrieve_gene_name_N_organism(uniprot_ref)
    myorganism=organism_tmp
    names='name:"'+"+".join(re.findall("(?<=name:).*?(?=\+OR)",myname_tmp))+'"'
    if re.search("(?<=gene:).*",myname_tmp):
        gene='gene:"'+re.findall("(?<=gene:).*",myname_tmp)[0]+'"'
    #myname=names+gene
        myname='('+names+"+OR+"+gene+')'
    else:
        myname='('+names+')'
    myquery0='http://www.uniprot.org/uniprot/?query=' + 'organism:"' + myorganism.replace(" ","+") + '"+AND+' + myname
    myquery = myquery0 + '&sort=score&columns=id&format=tab'
    f = urllib.urlopen(myquery)           
    myfile0 = f.read()
    myquery2 = myquery0 + '&sort=score&columns=id,entry,name,sequence&format=fasta'
    f = urllib.urlopen(myquery2)           
    myfile = f.read()  
    f.close()
    FASTA = myfile.split('\n>')[1:]
    uniprotFASTA=[]
    uniprotID=[]
    dic_fasta={}
    c=0
    uniprotFASTA=[]
    for i in FASTA:
#        if  "".join(i.split("\n")[1:]) in dic_fasta and c>300:
#            pass
        if c>1000:
            pass
        else:
            x=i.split("\n")
            seq="".join(x[1:])
            res=">"+x[0]+"\n"+seq
            uniprotFASTA.append(res)
            uniprotID.append(i.split("|")[1])
            dic_fasta["".join(i.split("\n")[1:])]=0
            c+=1
    write_file(folder+'UNIPROT-ID.txt',uniprotID)
    #ws.write_file(folder+'UNIPROT-FASTA.fasta',uniprotFASTA)
    f=open(folder+'UNIPROT-FASTA.fasta',"w")
    f.write(">refseq\n"+refseq+"\n")
    f.write("\n".join(uniprotFASTA))
    f.close()
    FASTA_headers = ''
    FASTA_headers = [str(((uniprotFASTA[i]).splitlines()[0]+'\n')) for i in range(len(uniprotFASTA))]
    FASTA_headers = str(FASTA_headers).replace(', ','').replace("'",'').strip('[]').replace('\\n','\n') 
    fo = open(folder+'UNIPROT-FASTA-header.txt','w')       
    fo.write(FASTA_headers)
    return (uniprotID, uniprotFASTA,organism_tmp,names.replace("+"," ")[6:-1])
 

def retrieve_all_sequences_fasta(folder,myname,myorganism,refseq):
    myquery0='http://www.uniprot.org/uniprot/?query=' + 'organism:"' + myorganism.replace(" ","+") + '"+AND+name:"' + myname.replace(" ","+")+'"'
    myquery = myquery0 + '&sort=score&columns=id&format=tab'
    print myquery0
    f = urllib.urlopen(myquery)           
    myfile0 = f.read()
    
    myquery2 = myquery0 + '&sort=score&columns=id,entry,name,sequence&format=fasta'
    f = urllib.urlopen(myquery2)           
    myfile = f.read()  
    f.close()
    FASTA = myfile.split('\n>')[1:]
    uniprotFASTA=[]
    uniprotID=[]
    dic_fasta={}
    c=0
    uniprotFASTA=[]
    for i in FASTA:
#        if  "".join(i.split("\n")[1:]) in dic_fasta and c>300:
#            pass
        if c>500:
            pass
        else:
            x=i.split("\n")
            seq="".join(x[1:])
            res=">"+x[0]+"\n"+seq
            uniprotFASTA.append(res)
            uniprotID.append(i.split("|")[1])
            dic_fasta["".join(i.split("\n")[1:])]=0
            c+=1
    write_file(folder+'UNIPROT-ID.txt',uniprotID)
    #ws.write_file(folder+'UNIPROT-FASTA.fasta',uniprotFASTA)
    f=open(folder+'UNIPROT-FASTA.fasta',"w")
    f.write(">refseq\n"+refseq+"\n")
    f.write("\n".join(uniprotFASTA))
    f.close()
    FASTA_headers = ''
    FASTA_headers = [str(((uniprotFASTA[i]).splitlines()[0]+'\n')) for i in range(len(uniprotFASTA))]
    FASTA_headers = str(FASTA_headers).replace(', ','').replace("'",'').strip('[]').replace('\\n','\n') 
    fo = open(folder+'UNIPROT-FASTA-header.txt','w')       
    fo.write(FASTA_headers)
    return (uniprotID, uniprotFASTA,myorganism)   
    
def retrieve_all_sequences(folder,uniprot_ref):

    (myname,myorganism) = retrieve_gene_name_N_organism(uniprot_ref)            
    ## UNIPROT IDs
    myquery = 'http://www.uniprot.org/uniprot/?query=' + 'organism:' + myorganism + '+AND+%28' + myname + '%29&sort=score&columns=id&format=tab'
    f = urllib.urlopen(myquery)           
    myfile0 = f.read()   
    uniprotID = myfile0.splitlines()[1:] # remove first row where it is only written "Entry"
    np.save(folder+'UNIPROT-ID',uniprotID)
    refidx = uniprotID.index(uniprot_ref) 
    

           
    ## FASTA sequence data
    myquery2 = 'http://www.uniprot.org/uniprot/?query=' + 'organism:' + myorganism + '+AND+%28' + myname + '%29&sort=score&columns=id,entry,name,sequence&format=fasta'
    f = urllib.urlopen(myquery2)           
    myfile = f.read()
    FASTA = myfile.split('>')[1:]
    
    uniprotFASTA=[]
    uniprotID=[]
    dic_fasta={}
    c=0
    for i in FASTA:
        if  "".join(i.split("\n")[1:]) in dic_fasta or c>300:
            pass
        else:
            uniprotFASTA.append(i)
            uniprotID.append(i.split("|")[1])
            dic_fasta["".join(i.split("\n")[1:])]=0
            c+=1    
    
    FASTA_headers = ''
    FASTA_headers = [str(((uniprotFASTA[i]).splitlines()[0]+'\n')) for i in range(len(uniprotFASTA))]
    FASTA_headers = str(FASTA_headers).replace(', ','').replace("'",'').strip('[]').replace('\\n','\n') 
    
    ## reference sequence
    fasta_ref = uniprotFASTA[refidx]
    refseq = fasta_ref[fasta_ref.index('\n'):].replace('\n','')

    ## Write files  
    write_file(folder+'UNIPROT-ID.txt',uniprotID)
    write_file(folder+'UNIPROT-FASTA.fasta',uniprotFASTA)
    fo = open(folder+'UNIPROT-FASTA-header.txt','w')       
    fo.write(FASTA_headers)
    fo.close()
    
    return (refidx,uniprotID,uniprotFASTA,refseq,myorganism)
    
    
def align_all_sequences(folder):

    ## perform CLUSTALw alignment
    if True: # not os.path.isfile(folder+'UNIPROT-FASTA.aln'):
        print '\n-- Performing clustalw alignment --\n'
        clustalw_cline = ClustalwCommandline(cmd=ep.clustalw_loc, infile=folder+'UNIPROT-FASTA.fasta', outorder='input')
        stdout, stderr = clustalw_cline()
    if re.search("Only 1 sequence",stdout):
        f=open(folder+'UNIPROT-FASTA.fasta',"r")
        read=f.readlines()
        list_seq=np.array([[i for i in read[1][:-1]]])
        list_id=[read[0]]
    else:
        clustalw_align = AlignIO.read(folder+'UNIPROT-FASTA.aln', 'clustal')
        print('..... Alignment:\n')    
        
        ## retrieve aligned sequences to identify variants among the virus isolates
        Nalign = len(clustalw_align) # = np.shape(clustalw_align)[0]
        L = np.shape(clustalw_align)[1]
        list_seq = np.zeros((Nalign,L),dtype=np.string0)
        
        list_id = []
        for record in range(Nalign): # loop over each isolate (row)
            list_id.append(clustalw_align[record].description)
            list_seq[record,:] = (list(clustalw_align[record].seq))   
        np.save(folder+'list_seq_uniprot',list_seq)
    
    return (list_seq,list_id)
    
def compute_position_specific_variability(folder,refidx,list_seq,in_FP_EmblID=1):

    refseq = list_seq[refidx,:].tolist()
    if len(refseq) == 1: # if it's a list of list, take the first list (sometimes happens when converting arrays to list and reciprocally)
        refseq = refseq[0]
    indi_ces = [i for i, x in enumerate(refseq) if x != '-'] # Find not '-' in the reference sequence ()

    # remove refidx=0 in list_seq if the reference sequence was not part of FP_EmblID (keep_only_subgenotypes_from_fishpathogens())
    if not in_FP_EmblID:
        list_seq = np.delete(list_seq,refidx,axis=0)
        print '\n..... !! The reference protein was removed from the dataset to compute within-speices position-specific variabilities ('+str(np.shape(list_seq))+') !!'

    BL = sm.blosum95

    L = np.shape(list_seq)[1]
    Nalign = np.shape(list_seq)[0]

    uAA = []   
    residue_conserved = [] # 1 if position consists of >95% likely mutations, 0 otherwise
    percent_conservation = [] # percent of conserved residues at each position
    av_score_conservation = []  # average pair-wise blosum score at each position (for all pairs of within-species sequences)
    list_variants_freq = []
    list_variants_YN = []
    
    # determine background frequency of the dataset        
    list_aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X','B','Z','-']
    pssm = np.zeros((len(list_aa),L),dtype=float)
    aa_dataset = list_seq.tolist()[0]
    [aa_dataset.extend(list_seq.tolist()[i]) for i in range(1,Nalign)]
    size_dataset = float(len(aa_dataset))
    bg = []
    for aa in list_aa:
        bg.append(float(aa_dataset.count(aa))/size_dataset)    
    
    for idxL in range(L):

        # determine frequency of the most common amino acid at each position
        uAA = list(set((list_seq[:,idxL]))) # unique amino acids  

        tmp = []
        for jdxL in range(len(uAA)):
            count_uAA = list(list_seq[:,idxL]).count(uAA[jdxL])
            if uAA[jdxL] != '-':
                tmp.append(count_uAA)
            idx_aa = list_aa.index(uAA[jdxL])         
            
            if count_uAA > 0:
                pssm[idx_aa,idxL] = (count_uAA/float(Nalign))/bg[idx_aa]
            else: # this avoids the case 0/0
                pssm[idx_aa,idxL] = 0
        count_und = list(list_seq[:,idxL]).count('-') # for the percentage, remove '-' from the total count
        list_variants_freq.append(int(100*max(tmp)/float(Nalign-count_und)))
        if int(100*max(tmp)/float(Nalign-count_und)) < 95:
            list_variants_YN.append(1) # i.e. there is a variant
        else:
            list_variants_YN.append(0)
        
        # determine score of each position using BLOSUM95 matrix    
        SEQ = str((list_seq[:,idxL]).tolist()).replace("'","").replace(",","").replace(' ','').strip('[]')    # positions were reference sequences has no '-' (useful for visualization in Protter)         
        score_tmp = []
        for jdxL2 in range(len(SEQ)): # pair-wise computation of BLOSUM score
            for kdx in range(jdxL2+1,len(SEQ),1): # ensures no "self-self" scores, e.g. seq1 with seq1
                try:
                    score_tmp.append(BL[SEQ[jdxL2],SEQ[kdx]])
                except: # if '-'
                    score_tmp.append(BL[SEQ[kdx],SEQ[jdxL2]])
        av_score_conservation.append(np.mean(score_tmp))
        neg_scores = [i for i, x in enumerate(score_tmp) if x < 0]
        non_gapped_pairs = [i for i,x in enumerate(score_tmp) if not np.isnan(x)]
        N_pairs = len(non_gapped_pairs)
        if N_pairs>0:
            percent_conservation.append(1-float(len(neg_scores))/N_pairs)           
        else:
            percent_conservation.append(float('nan'))
        if len(neg_scores) > 0.05*N_pairs: # more than 5% negative scores (remove - to calculate)
            residue_conserved.append(len(neg_scores))
        else:
            residue_conserved.append(0)        
    
    # pssm   
    print('\n..... PSSM:\n')
    pssm = np.log(pssm)
    (fig,ax) = uf.new_figure(size=(12,7))
    ax.xaxis.set_tick_params(labeltop='on',labelbottom='off')
    plt.xticks(range(0,L,50))
    ax.set_yticklabels(list_aa)
    plt.yticks(range(0,len(list_aa)))
    plt.imshow(pssm,interpolation='nearest',cmap=plt.cm.coolwarm,aspect='auto')
    plt.colorbar()
    plt.show()
    fig.savefig(folder+'PSSM.pdf', dpi=fig.dpi)
    
    # keep only those position with no gap ('-') in the "reference" sequence
    percent_conservation = [percent_conservation[i] for i in indi_ces] 
    av_score_conservation = [av_score_conservation[i] for i in indi_ces] 
    score = [residue_conserved[i] for i in indi_ces]
    list_variants_YN = [list_variants_YN[i] for i in indi_ces]

    indices_within_species_variability = uf.indices2string([i for i, x in enumerate(score) if x > 0])           
    #indices_within_species_conservation = uf.indices2string([i for i, x in enumerate(score) if x == 0])  

    fo = open(folder+'indices_within_species_variability.txt','w')        
    fo.write(str(indices_within_species_variability))
    fo.close()    
    
    return (percent_conservation,av_score_conservation,pssm)

def get_ws_output():
    
    ws_output = np.load('within_species_output.npz')
    refidx = int(ws_output['refidx'])
    ws_uniprotID = ws_output['ws_uniprotID'].tolist()
    ws_uniprotFASTA = ws_output['ws_uniprotFASTA'].tolist()
    myorganism = str(ws_output['myorganism'])
    list_ws_seq = ws_output['list_ws_seq']
    list_ws_id= ws_output['list_ws_id'].tolist()
    percent_ws_conservation = ws_output['percent_ws_conservation'].tolist()
    av_score_ws_conservation = ws_output['av_score_ws_conservation'].tolist()
    ws_pssm = ws_output['ws_pssm']
    
    return (refidx,ws_uniprotID,ws_uniprotFASTA,myorganism,list_ws_seq,list_ws_id,percent_ws_conservation,av_score_ws_conservation,ws_pssm)
    
def get_conserved_positions(list_seq):
    
    S = np.shape(list_seq)
    N = float(S[0])
    L = S[1]

    conserved_positions = []
    percent_conserved = []
    
    for i in range(L):
        if list_seq[0][i]!="-":
            seq = list(list_seq[:,i])
            nn = seq.count('-')
            if nn>0:
                [seq.remove('-') for j in range(nn)]
            u_aa = list(set(seq))
            conserved_positions.append(1) if len(u_aa)==1 else conserved_positions.append(0)
            perc = 0
            for aa in u_aa:
                perc_tmp = seq.count(aa)/(N-nn)
                if perc_tmp>perc:
                    perc = perc_tmp 
            percent_conserved.append(perc)
            
    return (conserved_positions,percent_conserved)
