"""
@author: aitana neves | heig-vd
@author: jibril mammeri | heig-vd
"""


from robobrowser import RoboBrowser
import os, time, re, subprocess
from shutil import copyfile
from bs4 import BeautifulSoup
from Bio import Entrez, SeqIO
Entrez.email="xavier.brochet@heig-vd.ch"
import sys
sys.path.insert(0, '/code/westbesel/static/prediction/scripts/')
sys.path.insert(0, '/code/westbesel/static/prediction/scripts/epitope_prediction/')
import thread
from threading import Thread
import numpy as np
import epitope_prediction as ep
import xlsxwriter
import directories
import within_species as ws
import between_species as bs
import protein_characterization as pc
import predict_epitopes as pe
import html2text
import smtplib
import csv
import urllib2
import email.message
import matplotlib
import pickle
from collections import OrderedDict
from math import *
matplotlib.use('Agg')

SVMTrip_dir='/code/westbesel/static/prediction/scripts/SVMTriP/SVMTriP/SVMTriP/'

repeat=True
cystein=True


def write_excel_file(folder,project_name,ml,refseq,hphil_cons_stretch,hphile_epitope,solvent_accessibility,jury_score,weighted_jury_score,max_repeat_length,pos_cons_stretch,L_stretch,idx_pareto,pareto_dominant):
    
    workbook = xlsxwriter.Workbook(folder+'LIST-final-epitopes.xls')
    format_heading = workbook.add_format({'bold': True})
    format_bold = workbook.add_format({'bold': True})
    format_bold.set_text_wrap()
    format_bold_color = workbook.add_format({'bold': True})
    format_bold_color.set_font_color('red')
    format_cell_color = workbook.add_format()
    format_bold.set_text_wrap()
    format_cell = workbook.add_format()
    format_bold.set_font_size(9)
    format_bold_color.set_font_size(9)
    format_cell.set_font_size(7)
    format_cell.set_align('center')
    format_cell.set_align('vcenter')       
    format_bold.set_align('center')
    format_bold.set_align('vcenter')  
    format_bold_color.set_align('center')
    format_bold_color.set_align('vcenter') 
    format_cell_color.set_font_size(7)
    format_cell_color.set_align('center')
    format_cell_color.set_align('vcenter')       
    format_cell_color.set_font_color('blue')
    format_cell_color.set_italic()
    worksheet = workbook.add_worksheet()
    worksheet.set_landscape()
    worksheet.set_column('A:A', 4)
    worksheet.set_column('B:B', 3)
    worksheet.set_column('C:D', 14)
    worksheet.set_column('E:G', 12)
    worksheet.set_column('H:L', 9)
    worksheet.merge_range('A1:C1',project_name,format_heading)
    worksheet.write(1,0,'Start',format_bold)
    worksheet.write(1,1,'End',format_bold)
    worksheet.write(1,2,'Fragment',format_bold)
    worksheet.write(1,3,'Longest\nconserved\nfragment',format_bold)
    worksheet.write(1,4,'Average\nhydrophilicity\nstretch',format_bold_color)
    worksheet.write(1,5,'Average\nhydrophilicity\nepitope',format_bold)
    worksheet.write(1,6,'Solvent\naccessibility',format_bold_color)
    worksheet.write(1,7,'jury-score',format_bold)
    worksheet.write(1,8,'weighted\njury-score',format_bold_color)
    worksheet.write(1,9,'Max\nrepeat\nlength',format_bold)
    worksheet.write(1,10,'Length\nstretch',format_bold)#_color)
    worksheet.write(1,11,'Ref. number',format_bold)
    worksheet.write(1,12,'Pareto\nscore',format_bold)
    #for i in range(len(ep_candidates)):
    for x,i in enumerate(idx_pareto):
        worksheet.write(x+2,0,ml[i][0]+1,format_cell)  # add +1 to start counting from 1 (more intuitive in xls file and similar to Blast output)
        worksheet.write(x+2,1,ml[i][-1]+1,format_cell) # add +1 to start counting from 1 (more intuitive in xls file and similar to Blast output)         
        worksheet.write(x+2,2,refseq[ml[i][0]:ml[i][-1]+1],format_cell)     
        worksheet.write(x+2,3,refseq[pos_cons_stretch[i][0]:pos_cons_stretch[i][-1]+1],format_cell_color)    
        worksheet.write(x+2,4,hphil_cons_stretch[i],format_cell)
        worksheet.write(x+2,5,hphile_epitope[i],format_cell)
        worksheet.write(x+2,6,solvent_accessibility[i],format_cell)
        worksheet.write(x+2,7,jury_score[i],format_cell)
        worksheet.write(x+2,8,weighted_jury_score[i],format_cell)
        worksheet.write(x+2,9,max_repeat_length[i],format_cell)
        worksheet.write(x+2,10,L_stretch[i],format_cell)
        worksheet.write(x+2,11,i,format_cell)
        worksheet.write(x+2,12,pareto_dominant[i],format_cell)
        format_cell.set_font_color('black')
    workbook.close()

def pipeline(fname,db,proc,refseq='',myorganism="",myname=""):
        """
        Select epitopes, read results from others prediction and compute it all
        """
#    try:
        global repeat
        start=0
        
        #%% DIRECTORIES
        cur_dir = directories.get_current_directory('__file__')
        folder_DATA = '/code/westbesel/static/prediction/output/'
        folder_RESULTS = '/code/westbesel/static/prediction/RESULTS_test/'
        
        #%% USER-SPECIFIED INPUT
        Hphob_scale = ep.scale # hydrophobicity scale
        
        sub_title = '_aaCS'
        
        #----------------------------------------------------------------------------------------------------------------------
        in_FP_EmblID = 1 # DO NOT MODIFY
        project_name = fname+sub_title # short name
        uniprot_ref = fname
        folder_species = folder_DATA+fname+'/'
        data_folder =  folder_DATA+fname+'/'
        #----------------------------------------------------------------------------------------------------------------------
        #%% PROJECT DIRECTORIES
        (mypathR,subf) = directories.set_directories(folder_RESULTS,project_name,ep.sf2)
        #%% LIST OF WITHIN-SPECIES SEQUENCES
        folder = subf['within-species']['path']
        print "within species"
        proc.step="within species"
        proc.save()
        if db == "uniprot" :
            (ws_uniprotID,ws_uniprotFASTA,myorganism,myname)    = ws.retrieve_all_sequences_uniprot(folder,uniprot_ref,myname,myorganism,refseq)
            proc.protein = myname
            proc.organism = myorganism.replace("+"," ")
            proc.save()
        else:
            (ws_uniprotID,ws_uniprotFASTA,myorganism)    = ws.retrieve_all_sequences_fasta(folder,myname,myorganism,refseq)
        (list_ws_seq,list_ws_id)                         = ws.align_all_sequences(folder)
        (conserved_pos,percent_conserved)                = ws.get_conserved_positions(list_ws_seq)
    
        folder = subf['between-species']['path']
        proc.step = "between species"
        proc.save()
        stderr = None
        blast_output = bs.retrieve_homologous_sequences(folder,refseq,myorganism)
        try:
            (list_homo_seq,indices_homo,stderr)           = bs.align_all_sequences(blast_output)
            (conserved_pos2,percent_conserved2)           = ws.get_conserved_positions(list_homo_seq)
            (_,_,indices_bs_variability)                  = bs.compute_position_specific_variability(folder,list_homo_seq,indices_homo)
        except:
            conserved_pos2=[1]*len(refseq)
            percent_conserved2=[1.0]*len(refseq)
    #    
        #%% PROTEIN CHARACTERIZATION
        proc.step="protein characterization"
        proc.save()
        if db=="uniprot" and proc.translation:
            (indices_TM,indices_SIGNAL,indices_PTM, indices_mut)     = pc.get_uniprot_info(uniprot_ref)
        else:
            (indices_TM,indices_SIGNAL,indices_PTM, indices_mut)     = ([],[],[],[])
        position_specific_hydrophilicity_score                       = pc.get_hydrophobicity(Hphob_scale,refseq)
        (indices_B,buried_score)                                     = pc.get_solvent_accessibility_JPred(data_folder)
        (indices_H,indices_E)                                        = pc.get_secondary_structure_JPred(data_folder)
        (indices_H_s,indices_E_s,indices_DR_s,indices_B_s)           = pc.get_scratch_predictions(data_folder)
        indices_H_consensus                                          = pc.get_consensus_predictions(indices_H,indices_H)
        indices_E_consensus                                          = pc.get_consensus_predictions(indices_E,indices_E)
        indices_B_consensus                                          = pc.get_consensus_predictions(indices_B,indices_B)
        #%% GET CANDIDATE EPITOPES DISCARD SOME EPITOPES TO CRITERIA & POST-PROCESS REMAINING ONES
        proc.step="predict epitopes"
        proc.save()
        sliding_windows                                              = pe.get_all_sliding_windows(refseq,indices_SIGNAL) # including shorter at terminal regions (fenetre de 20AA qui glisse sur la sequence)
        if proc.repeats:
            (ml,_)                                                   = pe.discard_epitopes_with_repeats(sliding_windows,refseq)
        else:
            ml                                                       = sliding_windows
        if proc.translation == True :
            ml                                                       = pe.discard_epitopes_with_PTM(ml,indices_PTM)
        ml                                                           = pe.discard_hydrophobic_epitopes(ml,position_specific_hydrophilicity_score)
        if proc.cystein:
            ml                                                       = pe.post_process_cysteins(ml,refseq,Hphob_scale) 
        
        #%% KEEP ONLY EPITOPES WITH (CONSERVED & HYDROPHILIC) STRETCHES OF AT LEAST ep.MIN_STRETCH_LENGTH aa
        start=1
        (ep_candidates,pos_cons_stretch,hphil_cons_stretch,L_stretch) = pe.keep_epitopes_with_conserved_hydrophilic_stretches(ml,percent_conserved,ep.THRES_CONSERVATION,position_specific_hydrophilicity_score)
        
        #%% VIEW EPITOPE PREDICTIONS FROM ALL METHODS
        
        folder = subf['epitope-predictions-from-online-methods']['path']
        (epitopes,ml_majority,score_majority)                        = pe.get_all_epitope_predictions(folder,data_folder,len(refseq),project_name)  
        weighted_score_res                                           = pe.get_consensus_epitopes(epitopes,len(refseq),score_majority,clustering_folder=folder_DATA,folder=folder)
        output                                                       = pe.view_epitope_predictions('all-prediction-methods',folder,project_name,epitopes,ml_majority,score_majority,indices_SIGNAL,indices_TM,indices_PTM,percent_conserved,indices_bs_variability,indices_H_consensus,indices_E_consensus,indices_DR_s,buried_score)
        
        #%% WEBLOGOS OF ALL CANDIDATE EPITOPES
        folder = subf['weblogos-epitopes']['path']
        #Modified by jibril
        #output                                                       = pe.write_epitope_ws_fasta(folder,list_ws_seq,ws_uniprotID,ep_candidates,genotype,extra_info)
        output                                                       = pe.write_epitope_ws_fasta(folder,list_ws_seq,ws_uniprotID,ep_candidates)
        output                                                       = pe.draw_epitope_ws_weblogos(folder,len(ep_candidates))
        #%% BLAST ALL CANDIDATE EPITOPES
        
        #folder = subf['blast-epitopes']['path']
        
        #max_homology_identity                                        = pe.blast_all_epitopes(folder,folder_species,ep_candidates,refseq) 
        
        #%% RANK EPITOPES ACCORDING TO A LIST OF CRITERIA
        proc.step="write results"
        proc.save()
        folder = subf['final-predictions']['path']
        buried_score[np.where(buried_score==1)] = 5
        buried_score[np.where(buried_score==2)] = 25
        buried_score[np.where(np.isnan(buried_score))] = 100
        max_repeat_length                                            = [max(pe.compute_repeat_sequences(refseq[mli[0]:mli[-1]+1])) for mli in ep_candidates]
        hphile_epitope                                               = [np.mean(position_specific_hydrophilicity_score[mli]) for mli in ep_candidates]
        solvent_accessibility                                        = [np.mean(buried_score[mli]) for mli in ep_candidates]
        jury_score                                                   = [np.mean(score_majority[mli]) for mli in ep_candidates]      
        jury_score                                                   = [(np.array(jury_score[i])+1)/2.0 for i in range(len(jury_score))] 
        weighted_jury_score                                          = [np.mean(weighted_score_res[mli]) for mli in ep_candidates]
        weighted_jury_score                                          = [(np.array(weighted_jury_score[i])+1)/2.0 for i in range(len(weighted_jury_score))] 
        data                                                         = zip(hphil_cons_stretch,L_stretch,solvent_accessibility,weighted_jury_score)
        data2                                                        = zip(hphil_cons_stretch,solvent_accessibility,weighted_jury_score)
        (idx_pareto,pareto_dominant)                                 = pe.compute_pareto_dominance(data2)
        output                                                       = write_excel_file(folder,project_name,ep_candidates,refseq,hphil_cons_stretch,hphile_epitope,solvent_accessibility,jury_score,weighted_jury_score,max_repeat_length,pos_cons_stretch,L_stretch,idx_pareto,pareto_dominant)
        dic_end=OrderedDict()
        list_pos=[]
        ep_candidates2=[[] for x in ep_candidates]
        #We save a pickle dictionnary that will be used to generate the results page on the web interface
        #Fill the dictionnary with the data for the web interface
        for i in xrange(0,len(idx_pareto)):
            ep_candidates2[idx_pareto[i]] = ep_candidates[i]
        for i in ep_candidates2:
            dic_end[str(i[0])+"-"+str(i[-1])] = False
            list_pos.append([])
        for i in dic_end:
            for j in list_pos:
                if dic_end[i] == False:
                    start = int(i.split("-")[0])
                    end = int(i.split("-")[1])+1
                    boo = True
                    for k in range(start,end):
                        if k in j:
                            boo = False
                    if boo == True:
                        dic_end[i] = True
                        for k in range(start,end):
                            j.append(k)
        while [] in list_pos:
            list_pos.remove([])        
        dic_methods = OrderedDict()
        #The key "westbesel" contains a list of AA for each position and if this AA is in at least an epitope or not
        #The key "epitopes" the positions of all epitopes
        dic_methods["WestBeSEL"] = [[x,False] for x in refseq]
        dic_methods["epitopes"] = list_pos
        for i in ep_candidates:
            for j in i:
                dic_methods["WestBeSEL"][j][1] = True
        for i in epitopes:
            dic_methods[i] = list(np.repeat("",len(refseq)))
            #print len(dic_methods[i])
            for j in epitopes[i]["ml"]:
                for k in j:
                    if k<len(dic_methods[i]):
                        if i == "SVMTriP":
                            dic_methods[i][k+1] = "X"
                        elif i == "AAP":
                            dic_methods[i][k-1] = "X"
                        else:
                            dic_methods[i][k] = "X"
        #Conservation and variability contains the percentage of conserved AA for this position
        dic_methods["Conservation"] = percent_conserved
        dic_methods["Variability"] = percent_conserved2
        #Secondary strcuture
        dic_methods["Secondary structure"] = list(np.repeat("c",len(refseq)))
        #Disorder
        dic_methods["Disorder"] = list(np.repeat("",len(refseq)))
        #Hydrophobicity according to the scale defined earlier
        dic_methods["Hydrophobicity"] = list(np.repeat("",len(refseq)))
        #Post translational modifications
        dic_methods["PTM"] = list(np.repeat("",len(refseq)))
        # Mutations
        dic_methods["Mutation"] = list(np.repeat("",len(refseq)))
        #We fill each of these dictionary entries with the data obtained in this pipeline
        for i in indices_mut:
            dic_methods["Mutation"][i] = "-"        
        for i in indices_PTM:
            dic_methods["PTM"][i] = "-"
        dic_methods["Solvent accessibility"] = []
        for i in xrange(0,len(position_specific_hydrophilicity_score)):
            if position_specific_hydrophilicity_score[i]<0:
                dic_methods["Hydrophobicity"][i] = "-"
            elif position_specific_hydrophilicity_score[i]>0:
                dic_methods["Hydrophobicity"][i] = "+"
            else:
                dic_methods["Hydrophobicity"][i] = "-"
        #print "4"
        for i in buried_score:
            dic_methods["Solvent accessibility"].append(float(i)/100.0)    
        if len(indices_H_consensus) > 0:
            for i in indices_H_consensus:
                if i != "":
                    dic_methods["Secondary structure"][int(i)-1] = "H"
        if len(indices_E_consensus) > 0:
            for i in indices_E_consensus:
                if i != "":
                    dic_methods["Secondary structure"][int(i)-1] = "E"
        if len(indices_DR_s) > 0:
            for i in indices_DR_s.split(","):
                if i != "":
                    dic_methods["Disorder"][int(i)-1] = "D"
        #We save the dictionary with pickle
        dic2 = open(data_folder+"dic2","w")
        pickle.dump(dic_methods,dic2)
        dic2.close()
        os.chdir(folder_RESULTS)    
        command = "tar -cf " + project_name+".tar " + project_name
        os.system(command)
        os.chdir(cur_dir)
        proc.over=True
        proc.save()
        #Send an email to the users when the job is over
        Email(proc.mail,proc.hash_id)
    
def scratch(filename):
    #Reproduce the scratch file with NetSurfP and Disembl
    f = open(filename,"r")
    fasta = f.readlines()
    f.close()
    res = []
    fasta = "".join(fasta[1:]).replace("\n","")
    print "DisEMBL"
    res.append(fasta)
    res.append(fasta)
    res.append(DisEmbl(res[0]))
    res.append(fasta)
    print "DisEMBL done"
    #res.append(pred[2])
    #Cleanup
    out = open(filename.replace(".fasta","/")+"SCRATCH.txt","w")
    out.write("\n".join(res).replace("\n\n","\n"))
    out.close()      


#%%
    
def jpred(filename):   
    #dir_out='/home/jibril/Desktop/westbesel-master/westbesel/Implementation/DATA/P19691/'
    f = open(filename,"r")
    header = f.readline()
    corpus = ""
    for i in f.readlines():
        if ">" not in i:
            corpus += i
    f.close()
    if len(corpus) < 800:
        print "case 1" 
        command = "perl ../scripts/jpredapi submit mode=single format=fasta seq="+\
        corpus.replace("\n","")+" email=name@domain.com" 
        x = command.split(" ")
        proc = subprocess.Popen(x,stdout=subprocess.PIPE)
        proc.wait()
        for line in iter(proc.stdout.readline,''):
            if re.search("perl.*silent",line):
                com = re.search("perl.*silent",line).group(0)
        command = com.split(" ")
        command[1] = '../scripts/jpredapi'
        command[5] = 'checkEvery=10'
        proc = subprocess.Popen(command,stdout=subprocess.PIPE)
        proc.wait()
        for line in iter(proc.stdout.readline,''):
            if re.search("[A-Z-a-z0-9_]*\/.+tar\.gz",line):
                path=re.search("[A-Z-a-z0-9_]*\/.+tar\.gz",line).group(0)
       
        com = "tar xzvf "+path+" -C "+path.split("/")[0]
        proc = subprocess.Popen(com.split(" "),stdout=subprocess.PIPE)
        proc.wait()
        res = path.split(".")[0]+".jnet"
        com = ['cp',res,filename.replace(".fasta","/")+'JPred.txt']
        proc = subprocess.Popen(com,stdout=subprocess.PIPE)
        proc.wait()
        com = ['rm','-rf',path.split("/")[0]+'/']
        proc = subprocess.Popen(com,stdout=subprocess.PIPE)
    else:
        name = filename.split(".")[0].split("/")[-1]
        x = ceil(len(corpus)/600.0)
        jnet = [""]*12
        for i in xrange(0,int(x)):
            f = open("jpred"+name+str(i)+".fasta","w")
            f.write(header)
            f.write(corpus[i*600:(i+1)*600])
            f.close()
            command = "perl ../scripts/jpredapi submit mode=single format=fasta file="+\
            "jpred"+name+str(i)+".fasta"+" email=toto@thing.com" 
            x = command.split(" ")
            proc = subprocess.Popen(x,stdout=subprocess.PIPE)
            proc.wait()
            for line in iter(proc.stdout.readline,''):
                if re.search("perl.*silent",line):
                        com = re.search("perl.*silent",line).group(0)
            command = com.split(" ")
            command[1] = '../scripts/jpredapi'
            command[5] = 'checkEvery=10'
            proc = subprocess.Popen(command,stdout=subprocess.PIPE)
            proc.wait()
            for line in iter(proc.stdout.readline,''):
                if re.search("[A-Z-a-z0-9_]*\/.+tar\.gz",line):
                    path = re.search("[A-Z-a-z0-9_]*\/.+tar\.gz",line).group(0)
            com = "tar xzvf "+path+" -C "+path.split("/")[0]
            proc = subprocess.Popen(com.split(" "),stdout=subprocess.PIPE)
            proc.wait()
            res=path.split(".")[0]+".jnet"
            #print res
            f = open(res,"r")
            lines = f.readlines()
            for j in xrange(1,len(lines)):
                if i == 0:
                    jnet[j] = lines[j].replace("\n","")
                else:
                    jnet[j] = jnet[j]+lines[j].split(":")[1].replace("\n","")
            com = ['rm','-rf',path.split("/")[0]+'/']
            proc = subprocess.Popen(com,stdout=subprocess.PIPE)
        out = open(filename.replace(".fasta","/")+'JPred.txt',"w")
        out.write("\n".join(jnet))
        out.close()
        
#%%
def get_fasta_from_AN(AN):
    handle = Entrez.efetch(id=AN,db='protein',rettype='fasta',retmode='text')
    record = SeqIO.read(handle, "fasta")
    fasta = ">"+record.description+'\n'+str(record.seq)
    f = open(record.id+".fasta","w")
    f.write(fasta)
    f.close()
    return record.id+".fasta"
#%%
def get_fasta_from_uniprot_ID(ID):
    command = "wget "+"http://www.uniprot.org/uniprot/"+ID+".fasta"
    os.system(command)
    return ID+".fasta"
#%%
def create_directory(name):
    command = "mkdir "+name
    os.system(command)
#%%
def Epitope_prediction(fasta_file, name):
    fasta=""
    f = open(fasta_file,"r")
    for i in f:
        fasta += i
    f.close()
    name = name+"/"
    t = []
    t.append(Thread(target=BCPPred_AAP, args=(fasta,"bcpred",name)))
    t.append(Thread(target=BCPPred_AAP, args=(fasta,"aap",name)))
    t.append(Thread(target=BepiPred, args=(fasta_file,name)))
    t.append(Thread(target=ABCPred, args=(fasta_file,name)))
    t.append(Thread(target=scratch, args=(fasta_file,)))
    for i in t:
        i.start()
    for i in t:
        i.join()
    print "Jpred Start"
    jpred(fasta_file)
    print "Jpred end"
    SVMTriP(fasta_file,name)
    os.system("rm "+fasta_file)
#%%
def BCPPred_AAP(fasta,method,name,specificity='75'):
    fasta='\r\n'.join(fasta.split("\n")[1:])
    for length in ['12','20','16']:
        try :
            browser = RoboBrowser()
            login_url = 'http://ailab.ist.psu.edu/bcpred/predict.html'
            browser.open(login_url)
            form = browser.get_forms()
            form[1]["sequence"].value = fasta
            form[1]['pmethod'].options = ["bcpred","aap"]
            form[1]["pmethod"].value = method
            form[1]["length"].value = length
            form[1]["specificity"].value = specificity
            browser.submit_form(form[1])
            x = browser.response.content
            y = re.search('<table(.+\n)+<\/P>',x)
            z = y.group(0)
            soup = BeautifulSoup(z)
            x=soup.get_text().split('\n')[2:]
            if method=="aap":
                out=open(name+"AAP_"+length+"aa_NO.txt","w")
            else:
                out=open(name+"BCPred_"+length+"aa_NO.txt","w")
            for i in x:
                if i=="":
                    break
                l=i.replace(" ","")
                line=""
                next=0
                for j in l:
                    if next==0 and re.search('[A-Z]',j):
                        line+="\t"
                        next=1
                    if next==1 and re.search('[0-9]',j):
                        line+="\t"
                        next=0
                    line+=j
                out.write(line+"\n")
            out.close()
            browser = RoboBrowser()
            login_url = 'http://ailab.ist.psu.edu/bcpred/predict.html'
            browser.open(login_url)
            form = browser.get_forms()
            form[1]["sequence"].value=fasta
            form[1]['pmethod'].options=["bcpred","aap"]
            form[1]["pmethod"].value=method
            form[1]["length"].value=length
            form[1]["specificity"].value=specificity
            form[1]["overlap"].value=[]
            browser.submit_form(form[1])
            x = browser.response.content
            y = re.search('<table(.+\n)+<\/TABLE>',x)
            z = y.group(0)
            soup = BeautifulSoup(z)
            x = soup.get_text().split('\n')[2:]
            if method == "aap":
                out=open(name+"AAP_"+length+"aa.txt","w")
            else:
                out=open(name+"BCPred_"+length+"aa.txt","w")
            for i in x:
                if i == "":
                    break
                l = i.replace(" ","")
                line = ""
                next = 0
                for j in l:
                    if next == 0 and re.search('[A-Z]',j):
                        line += "\t"
                        next = 1
                    if next == 1 and re.search('[0-9]',j):
                        line+="\t"
                        next = 0
                    line+=j
                out.write(line+"\n")
            out.close()
        except:
            if method == "aap":
                out=open(name+"AAP_"+length+"aa.txt","w")
            else:
                out=open(name+"BCPred_"+length+"aa.txt","w")
            out.close()
    print method, "done" 

def DisEmbl(fasta):
    browser = RoboBrowser()
    url="http://dis.embl.de/"
    browser.open(url)
    form = browser.get_forms()
    form[0]["sequence_string"] = fasta
    form[0]["doApplet"] = []
    browser.submit_form(form[0])
    time.sleep(2)
    x = browser.response.content
    show = 0
    res=""
    for i in x.split("\n"):
        if show > 0:
            show += 1
        if re.search("REM465",i):
            show += 1
        if re.search("</p>",i):
            show = 0
        if show > 2:
            res = i
    res = re.sub("<.*?>","",res)
    res = res.replace(" ","")
    seq = ""
    for i in res:
        if i.isupper():
            seq += "D"
        else:
            seq += "O"
    if seq == "":
        seq = "O"*len(fasta)
    return seq


#%%
def BepiPred_one(fasta,name,threshold=0.35):
    try:
        fasta = fasta.replace("\n",'\r\n')
        browser = RoboBrowser()
        login_url = 'http://www.cbs.dtu.dk/services/BepiPred-1.0/'
        browser.open(login_url)
        form = browser.get_forms()
        form[0]["SEQPASTE"].value = fasta
        form[0]["threshold"].value = threshold
        browser.submit_form(form[0])
        #TO IMPROVE
        time.sleep(10)
        x = browser.response.content
        y = x.split("\n")
        out = open(name+"BepiPred.txt","w")
        for i in y:
            if re.search('.*(bepipred-1.0b epitope).*',i):
                out.write(i+"\n")
        out.close()
        print "Bepipred done"
    except:
        #s.kill()
        print "Job error "+name
        out = open(name+"BepiPred.txt","w")
        out.close()


def BepiPred(fasta_file, name):
    try:
        browser = RoboBrowser()
        login_url = 'http://www.cbs.dtu.dk/services/BepiPred-2.0/'
        browser.open(login_url)
        form = browser.get_forms()
        form[0]["uploadfile"].value = open(fasta_file,"r")
        browser.submit_form(form[0])
        #TO IMPROVE
        time.sleep(3)
        x="Send me email when job finishes"
        login_url = browser.url
        browser.open(login_url)
        while (x in browser.response.content):
            login_url = browser.url
            browser.open(login_url)
            time.sleep(3)
        login_url = browser.url
        browser.open(login_url)
        ID=login_url.split("=")[1].split("&")[0]
        url='http://www.cbs.dtu.dk/services/BepiPred-2.0/tmp/'+ID+'/summary_'+ID+'.csv'
        response = urllib2.urlopen(url)
        cr = csv.reader(response)
        out = open(name+"BepiPred.txt","w")
        for row in cr:
            out.write("\t".join(row)+"\n")
        out.close()
        print "Bepipred done"
    except Exception as e:
        #s.kill()
        print "Job error "+name
        print e
        out = open(name+"BepiPred.txt","w")
        out.close()

def SVMTriP(fasta_file,name):
    try:
        directory = os.getcwd()
        dest = SVMTrip_dir
        os.system("cp "+fasta_file+' '+dest)
        os.chdir(dest)
        fasta = ">000\n"    
        f = open(fasta_file.split("/")[-1],"r")
        for i in f:
            if not re.search(">",i):
                fasta += i
        f.close()
        f = open(fasta_file.split("/")[-1],"w")
        f.write(fasta)
        f.close()
        print "Job start "+name
        #s=subprocess.Popen(["./svm_classify_server", "model_16aa.dat", "8001"])
        p = subprocess.Popen(["perl", "SVMTriP.pl", fasta_file.split("/")[-1], "16"])
        p.wait()
        #s.kill()
        results=[]
        f = open(fasta_file.split("/")[-1]+"_16aa_optimal_predictions.txt","r")
        for i in f:
             results.append(i)
        f.close()
        #subprocess.Popen(["rm ",fasta_file.split("/")[-1]+"*"])
        os.chdir(directory)
        print "Job finished "+name
        out=open(name+"SVMTriP.txt","w")
        for i in xrange(0,len(results)):
            j = results[i].replace(" ","")
            pos,seq,score = j.split(",")
            score = score.split('\t')[0]
            pos1,pos2 = pos.split("-")
            line = "\t".join([str(i+1),pos1,pos2,seq,score+"\n"])
            out.write(line)
        out.close()
        print "SVMTrip done"
    except Exception as e:
        #s.kill()
        print "Job error "+name
        print e
        os.chdir(directory)
        out = open(name+"SVMTriP.txt","w")
        out.close()
        

#%%
def ABCPred(fasta_file,name):
    f = open(fasta_file,"r")
    lines = f.readlines()
    f.close()
    fasta = "\n".join(lines[1:])
    browser = RoboBrowser()
    login_url = 'http://crdd.osdd.net/raghava/abcpred/ABC_submission.html'
    browser.open(login_url)
    form = browser.get_forms()
    form[0]["SEQ"].value = fasta
    form[0]["Threshold"].value = '0.35'
    browser.submit_form(form[0])
    #TO IMPROVE
    time.sleep(10)
    x = browser.response.content
    res = re.search('<table BORDER=7 CELLSPACING=3 CELLPADDING=2 cols=2 WIDTH="75% bgcolor="#e5e5e5><TR><TD WIDTH="10%">Rank.*table>',x)
    x = res.group(0)
    res = html2text.html2text(x)
    res2 = res.replace("\n","").split("|")
    x=[]
    for i in xrange(6,201,5):
        x.append('\t'.join(res2[i:i+4]).replace(" ","").replace("---",""))
    out = open(name+"ABCPred.txt","w")
    out.write("\n".join(x))
    out.close()
    print "ABCpred done"
#%%    

def Email(mail, hash_id):
    FROM = 'WeStBESel_server@heig-vd.ch'
    TO = [mail]
    m = email.message.Message()
    m['From'] = FROM 
    m['To'] = mail
    m['Subject'] = "Epitope prediction results"
    m.set_payload("Your results are available at http://WeStBESel.iict.ch/results/"+hash_id);
    message = m.as_string()
    server = smtplib.SMTP("smtp.heig-vd.ch")
    server.sendmail(FROM, TO, message)
    server.quit()    

    
def prediction():
    global q
    while True:
        print "Len1",q.qsize()
        while not q.empty():
            print "Len2", q.qsize()
            args = q.get()
            print "Len3", q.qsize()
            try:
                print "Go"
                ID = args[0]
                database = args[1]
                proc = args[2]
                organism = args[3]
                start = "/code/westbesel/"
                output = "/code/westbesel/static/prediction/output/"
                os.chdir(output)
                if database == "uniprot":
                    filename = get_fasta_from_uniprot_ID(ID)
                    f = open(filename,"r")
                    seq = "".join(f.readlines()[1:]).replace("\n","")
                    f.close()
                    create_directory(ID)
                    print "START"
                    Epitope_prediction(filename, ID)
                    print "END"
                    #subprocess.Popen(["rm ",filename])
                    os.chdir(start)
                    pipeline(ID,database,proc,seq,organism)
                elif database == "ncbi":
                    filename = get_fasta_from_AN(ID)
                    create_directory(ID)
                    print "START"
                    Epitope_prediction(filename, ID)
                    print "END"
                    f = open(filename,"r")
                    seq = f.readlines()
                    f.close()
                    subprocess.Popen(["rm ",filename])
                    os.chdir(start)
                    pipeline(ID,database,proc,seq,organism)
            #    elif database=="ncbi":
            #        filename=get_fasta_from_AN(ID)
                else:
                    print "Available database are 'uniprot' and 'ncbi'"
            except Exception as e:
                print e
                proc = args[2]
                proc.over = True
                proc.fail = str(e)
                print str(e)
                proc.save()
            print "Len4", q.qsize()
            q.task_done()
            print "Len5", q.qsize()
        time.sleep(30)
        print "Queue is empty, thread still running, no issue to report"
    
#%%    
        
if __name__ == '__main__':
    prediction(sys.argv[1],sys.argv[2])
#fasta_file="/home/jibril/Desktop/Epitope_predictors/P19691.fasta"
#fasta="MTSALRETFTGLRDIKGGVLEDPETEYRPGTITLPLFFSKADFDLEMIKRAVSHVGGEGTRRALGLLCAFVIAETVPSGRGTVAELLEALGFLLESLETGAPLEVTFADPNNKLAETIVKENVLEVVTGLLFTCALLTKYDVDKMATYCQNKLERLATSQGIGELVNFNANRGVLARIGAVLRPGQKLTKAIYGIILINLSDPATAARAKALCAMRLSGTGMTMVGLFNQAAKNLGALPADLLEDLCMKSVVESARRIVRLMRIVAEAPGVAAKYGVMMSRMLGVGYFKAYGINENARITCILMNINDRYDDGTSGGLTGLKVSDPFRKLAREIARLLVLKYDGDGSTGEGASDLIRRAEMASRGPDMGEEEEEDEEDDDSSEPGDSDSFL"
