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


