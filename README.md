
WeStBESel is a tool to help selecting the most relevant B-cell epitopes according to your needs (vaccines, diagnostic test...). We have used a mixed approach by combining physico-chemical property-based approaches, sliding windows algorithm combined with a panel of 5 of the best B-cell epitopes prediction tools based on machine learning methodologies.

You can either submit a FASTA sequence with the protein and organism name or you can use an Uniprot ID, in which case the sequence features will be automatically retrieved. You will receive an email when your results are available.

The WeStBESel output is divided into three sections, at the top of the page is the summary of the parameters used to analyze the target protein sequence, below a graphical output and following by is a summary of the results in tabular form. The results of the analysis are distributed in three tabs, you can use the navigation bar to flip through the various output pages, either the summary results (default page) or the multiple alignment results (intra-species or inter-species). 

############################<br>
<p>        INSTALLATION</p>
############################<br>

-Install docker (https://docs.docker.com/install/)<br>
-Unzip archive<br>
-Build docker image<br>
-Launch<br>

<br><br>
############################<br>
          DOCKER<br>
############################<br>
<br>
Build image:<br>
docker build -t ci4cb/westbesel . <br>
<br>
Launch docker:<br>
docker run -d -p 80:80 --env DJANGO_PRODUCTION=false ci4cb/westbesel<br>
<br>


