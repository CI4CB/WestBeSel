############################
        INSTALLATION
############################

-Install docker (https://docs.docker.com/install/)
-Unzip archive
-Build docker image
-Launch


############################
          DOCKER
############################

Build image:
docker build -t ci4cb/westbesel . 

Launch docker:
docker run -d -p 80:80 --env DJANGO_PRODUCTION=false ci4cb/westbesel

See active container :
docker ps

----OUTPUT----
CONTAINER ID        IMAGE               COMMAND                  CREATED             STATUS              PORTS                                              NAMES
bdc80459a5d2        jibril/westbesel    "/usr/bin/supervisord"   About an hour ago   Up About an hour    2525/tcp, 3306/tcp, 0.0.0.0:80->80/tcp, 8000/tcp   clever_bhabha

Stop container :
docker stop [CONTAINER ID] 

----EXAMPLE----
docker stop bdc80459a5d2
docker stop bdc8

See logs of active container:
docker ps -q | awk '{print $1}' | xargs docker logs -f --tail 20;

############################
       DJANGO TEMPLATE
############################


Template directory path : 
westbesel-master/westbesel/prediction/templates/prediction/

base.html : Header, footer, design
index.html : Submission form
results.html : Results page

Media path :
westbesel-master/westbesel/static/prediction/media



