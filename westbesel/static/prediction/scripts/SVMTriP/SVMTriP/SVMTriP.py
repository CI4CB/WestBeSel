import subprocess
import os

os.chdir("/code/westbesel/static/prediction/scripts/SVMTriP/SVMTriP/SVMTriP/")
subprocess.Popen(["./svm_classify_server", "model_16aa.dat", "8001"])
p=subprocess.Popen(["perl", "SVMTriP.pl", "P19691.fasta", "16"])
p.wait()
f=open("P19691.fasta_16aa_optimal_predictions.txt","r")
for i in f:
    print i
f.close()
