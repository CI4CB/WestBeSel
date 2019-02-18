# -*- coding: utf-8 -*-

"""
@author: jibril mammeri | heig-vd
"""

from __future__ import unicode_literals
from django.shortcuts import render, redirect
from django.templatetags.static import static
from django.http import HttpResponse
from django.contrib.auth import login, authenticate
from django.contrib.auth.forms import UserCreationForm
from django import forms
from prediction.models import *
from django.contrib.auth.models import User
import requests
import pandas as pd
import epitope_tools as ept
from threading import Thread
import hashlib
from Queue import Queue
global cowan_whittaker_Hphob_scale
import time, os
import pickle
from datetime import timedelta
from django.contrib.auth import logout
from django.http import Http404



queue=[]

def job_queue():
    """
    Queue to handle the jobs
    """
    global queue
    while True:
        if len(queue)>0:
            print queue
            args=queue[0]
            print args
            try:
                print "Go"
                ID=args[0]+"_"+str(args[2].id)
                #ID=args[0]                
                database=args[1]
                proc=args[2]
                ept.ep.THRES_CONSERVATION=proc.threshold
                ept.ep.MIN_STRETCH_LENGTH=proc.stretch
                ept.ep.min_epitope_length=proc.min_ep_length
                ept.ep.max_epitope_length=proc.max_ep_length
                ept.ep.TERMINAL_REGION=proc.max_ter_length
                ept.ep.MIN_TERMINAL_REGION=proc.min_ter_length
                organism=args[3]
                protein=args[4]
                fasta=args[5]
                start="/code/westbesel/"
                output="/code/westbesel/static/prediction/output/"
                os.chdir(output)
                if database=="uniprot":
                    filename=ept.get_fasta_from_uniprot_ID(ID)
                    f=open(filename,"r")
                    seq="".join(f.readlines()[1:]).replace("\n","")
                    f.close()
                    ept.create_directory(ID)
                    print "START"
                    proc.step="prediction methods"
                    proc.save()
                    ept.Epitope_prediction(filename, ID)
                    print "END"
                    os.chdir(start)
                    proc.refseq=seq
                    proc.step="pipeline start"
                    proc.save()
                    ept.pipeline(ID,database,proc,seq,organism)
                else:
                    ept.create_directory(ID)
                    fasta=fasta.replace("\r","")
                    f=open(ID+".fasta","w")
                    f.write(fasta)
                    f.close()
                    filename=ID+".fasta"
                    f=open(filename,"r")
                    seq="".join(f.readlines()[1:]).replace("\n","")
                    f.close()
                    print "START"
                    proc.step="prediction methods"
                    ept.Epitope_prediction(filename, ID)
                    print "END"
                    os.chdir(start)
                    proc.refseq=seq
                    proc.step="pipeline start"
                    proc.save()
                    ept.pipeline(ID,database,proc,seq,organism,protein)
                    #####       REPRENDRE ICI                     ####
                    
            except Exception as e:
                print e
                proc=args[2]
                proc.failures+=1
                print str(e)
                proc.save()
                if proc.failures >=3:
                    proc.over=True
                    proc.fail=str(e)
                else : 
                    if database == "uniprot":
                        queue.append((args[0],database,proc, "","",""))
                    else:
                        queue.append((args[0],database,proc, organism, protein, fasta))
            del queue[0]
        else:
            pass
            
worker=Thread(target=job_queue)
worker.setDaemon(True)
worker.start()

def IndexView(request):
    """
    View for the index page, form reading
    """
    global queue
    if request.POST:
        if len(request.POST.get('uniprot')) > 0 :
            name=request.POST.get('uniprot')
            db="uniprot"
            myorganism=""
            url = 'https://www.uniprot.org/uniprot/'+name
            r = requests.get(url)
            if not r.ok:
                return render(request, 'prediction/index.html',{"message":"Uniprot ID is not valid"})
        elif len(request.POST.get('protein')) > 0 and len(request.POST.get('fasta')) > 0 and len(request.POST.get('organism')) > 0:
            myprotein=request.POST.get('protein')
            myorganism=request.POST.get('organism')
            seq=request.POST.get('fasta')
            #if len(myorganism)< 5: 
            #    return render(request, 'prediction/index.html',{"message":"Please be as specific as possible for the organism name (example : Infectious hematopoietic necrosis virus (strain Round Butte))"})
            if seq[0] != ">":
                return render(request, 'prediction/index.html',{"message":"Fasta sequence first line must begin with a '>' "})
            for i in "".join(seq.split("\n")[1:]):
                if i not in "ACDEFGHIKLMNPQRSTVWY\r":
                    print "problem : "+i
                    return render(request, 'prediction/index.html',{"message":"Your sequence contains something other than amino acids "})
            name=myprotein.replace(" ",'-')
            db="fasta"
        scale=request.POST.get('hydrophob')
        #Settings reading
        if scale=="cowan":
            ept.ep.scale=ept.ep.cowan_whittaker_Hphob_scale
        else:
            ept.ep.scale=ept.ep.roseman_Hphob_scale
        email=request.POST.get('mail')
        proc=Process(name=name,mail=email,over=False,cystein=ept.cystein,repeats=ept.repeat )
        proc.step="start"
        proc.save()
        proc.hash_id=hashlib.md5(str(proc.id)).hexdigest()
        if request.user.is_authenticated():
            proc.user=request.user
            proc.exp_date=proc.date+timedelta(days=15)
        else:
            proc.exp_date=proc.date+timedelta(days=7)
        proc.save()
        proc.min_ep_length= int(request.POST.get('min_l'))
        proc.max_ep_length= int(request.POST.get('max_l'))
        proc.min_ter_length=int(request.POST.get('min_terminal_l'))
        proc.max_ter_length=int(request.POST.get('terminal_l'))
        proc.threshold= float(request.POST.get('threshold'))/100.0
        proc.stretch=int(request.POST.get('stretch_l'))
        if request.POST.get('cystein')=="yes":
            proc.cystein=True
        else:
             proc.cystein=False            
        if request.POST.get('repeat')=="yes":
            proc.repeat=True
        else:
            proc.repeat=False
        if request.POST.get('translation')=="yes":
            proc.translation=True
        else:
            proc.translation=False
        proc.save()
        #Case with an Uniprot ID
        if db=='uniprot':
            print "PUT"
            proc.uniprot = name
            proc.save()
            #q.put((name,db,proc, myorganism))
            queue.append((name,db,proc, "","",""))
            return redirect('results', proc_id=proc.hash_id)
            #Case with a FASTA sequence
        else:
            print "PUT"
            proc.organism = myorganism.replace("+"," ")
            proc.protein = myprotein
            proc.refseq = seq
            proc.save()
            queue.append((name,db,proc, myorganism, myprotein, seq))
            return redirect('results', proc_id = proc.hash_id)
    return render(request, 'prediction/index.html')

class SignUpForm(UserCreationForm):
    first_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
    last_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
    email = forms.EmailField(max_length=254, help_text='Required. Inform a valid email address.')

    class Meta:
        model = User
        fields = ('username', 'first_name', 'last_name', 'email', 'password1', 'password2', )


def UserView(request):
    """
    User page view
    """
    if request.user.is_authenticated():
        if request.method == 'POST':
            form=request.POST
            list_del=form.getlist("delete")
            for i in list_del:
                proc = Process.objects.get(id=i)
                name=proc.name+'_'+str(i)
                path1='/code/westbesel/static/prediction/output/'+name
                path2='/code/westbesel/static/prediction/RESULTS_test/'+name+"_aaCS"
                com='rm -rf '
                os.system(com+path1)
                os.system(com+path2)
                proc.delete()
            return redirect('user')
        username = request.user.username
        current_user = User.objects.get(username=username)
        proc=Process.objects.filter(user=current_user)
        return render(request, 'prediction/user.html',{"proc_list":proc})
    else:
        if request.method == 'POST':
            if 'login' in request.POST:
                form=request.POST
                username = form.get('login')
                raw_password = form.get('password')
                user = User.objects.get(username=username)
                success = user.check_password(raw_password)
                if success: 
                    user = authenticate(username=username, password=raw_password)
                    login(request, user)
                    return redirect('user')
                else:
                    return redirect('user')
            else:
                form = SignUpForm(request.POST)
                if form.is_valid():
                    form.save()
                    username = form.cleaned_data.get('username')
                    raw_password = form.cleaned_data.get('password1')
                    user = authenticate(username=username, password=raw_password)
                    login(request, user)
                    return redirect('user')
    form = SignUpForm()
    return render(request, 'prediction/user.html', {'form': form})

def HelpView(request):
    return render(request, 'prediction/help.html')


def excel(path_to_file,proc):
    """
    Extract data from an excel file to make a Python dictionnary
    """
    refseq=proc.refseq
    xl = pd.ExcelFile('/code/westbesel/static/prediction/RESULTS_test/'+proc.name+'_'+str(proc.id)+'_aaCS/final-predictions/LIST-final-epitopes.xls')
    df = xl.parse("Sheet1",header=1)
    dic={}
    for i in xrange(0,len(refseq)):
        dic[i+1]=[refseq[i],False]
    for i in xrange(0,len(df.Start)):
        for y in xrange(df.Start[i],df.End[i]+1):
            if dic[y][1]==False:
                dic[y][1]=True
    return dic

def logout_view(request):
    """
    Logout
    """
    logout(request)
    return redirect('index')            
            
def ResultView(request, proc_id=1):
    """
    View for the Results page
    """
    try:
        proc=Process.objects.get(hash_id=proc_id)
    except Process.DoesNotExist:
        raise Http404("Page does not exist")
    if proc.over and proc.fail==None:
        path_to_file = '/code/westbesel/static/prediction/RESULTS_test/'+proc.name+'_'+str(proc.id)+'_aaCS/final-predictions/'
        path_to_file2 = '/code/westbesel/static/prediction/output/'+proc.name+'_'+str(proc.id)+'/'
        xl = pd.ExcelFile(path_to_file+"LIST-final-epitopes.xls")
        df1 = xl.parse('Sheet1')
        #dic1=pickle.load(open(path_to_file+"dic1"))
        dic2=pickle.load(open(path_to_file2 + "dic2"))
        #print len(df1)
        return render(request, 'prediction/result.html',{"dic2":dic2,"Test":proc,"Sequence":excel(path_to_file+"LIST-final-epitopes.xls",proc),"path":"prediction/RESULTS_test//"+proc.name+'_'+str(proc.id)+"_aaCS/2_epitope-predictions-from-online-methods/FIGURE-all-prediction-methods.png","excel":df1[:101],"range":range(len(df1[:101]))})
    return render(request, 'prediction/result.html',{"Test":proc})

def Download(request, proc_id=1):
    """
    Download archive view
    """
    proc=Process.objects.get(hash_id=proc_id)
    path_to_file = '/code/westbesel/static/prediction/RESULTS_test/'+proc.name+'_'+str(proc.id)+'_aaCS.tar'
    f = open(path_to_file, 'r')
    response = HttpResponse(f, content_type='application/zip')
    response['Content-Disposition'] = 'attachment; filename=results_'+proc.name+'_'+str(proc.id)+'.tar'
    return response

