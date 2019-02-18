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
See active container :<br>
docker ps<br>
<br>
----OUTPUT----<br>
CONTAINER ID        IMAGE               COMMAND                  CREATED             STATUS              PORTS                                              NAMES<br>
bdc80459a5d2        jibril/westbesel    "/usr/bin/supervisord"   About an hour ago   Up About an hour    2525/tcp, 3306/tcp, 0.0.0.0:80->80/tcp, 8000/tcp   clever_bhabha<br>
<br>
Stop container :<br>
docker stop [CONTAINER ID] <br>
<br>
----EXAMPLE----<br>
docker stop bdc80459a5d2<br>
docker stop bdc8<br>
<br>
See logs of active container:<br>
docker ps -q | awk '{print $1}' | xargs docker logs -f --tail 20;<br>
<br>
############################<br>
       DJANGO TEMPLATE<br>
############################<br>
<br><br>

Template directory path : <br>
westbesel-master/westbesel/prediction/templates/prediction/<br>
<br>
base.html : Header, footer, design<br>
index.html : Submission form<br>
results.html : Results page<br>
<br>
Media path :<br>
westbesel-master/westbesel/static/prediction/media<br>



