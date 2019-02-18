FROM ubuntu:14.04

# Enable production settings by default; for development, this can be set to 
# `false` in `docker run --env`
ENV DJANGO_PRODUCTION=true

# Set terminal to be noninteractive
ENV DEBIAN_FRONTEND noninteractive

# Enable MySQL root user creation without interactive input
RUN echo 'mysql-server mysql-server/root_password password devrootpass' | debconf-set-selections
RUN echo 'mysql-server mysql-server/root_password_again password devrootpass' | debconf-set-selections

# Install packages
RUN apt-get update && apt-get install -y \
    git \
    libmysqlclient-dev \
    mysql-server \
    nginx \
    build-essential \
    make \
    python-dev \
    python-lxml \
    python-mysqldb \
    python-setuptools \
    clustalw \
    ghostscript \
    libstdc++6 \
    ncbi-blast+ \
    blast2 \
    libwww-perl \
    supervisor \
    wget \
    python-tk \
    gcc-4.9 \
    vim

RUN apt-get upgrade libstdc++6 -y

ENV DISPLAY :0
RUN easy_install pip

# Handle urllib3 InsecurePlatformWarning
RUN apt-get install -y libffi-dev libssl-dev libpython2.7-dev nano
RUN pip install urllib3[security] requests[security] ndg-httpsclient pyasn1


# Configure Django project
ADD . /code
ADD requirements.txt /code/requirements.txt
ADD libstdc++* /usr/lib/x86_64-linux-gnu/
RUN mkdir /djangomedia
RUN mkdir /static
RUN mkdir /gspr
RUN mkdir /logs
RUN mkdir /logs/nginx
RUN mkdir /logs/gunicorn
RUN mkdir /logs/supervisord
WORKDIR /code
RUN pip install -r /code/requirements.txt
RUN easy_install weblogo
#RUN cp libstdc++* /usr/lib/x86_64-linux-gnu/
WORKDIR /code
#RUN perl /code/westbesel/static/prediction/scripts/SCRATCH-1D_1.1/install.pl
#RUN chmod 7777 -R /code/westbesel/prediction/scripts/SCRATCH-1D_1.1/tmp/
#RUN umask 000 -R /code/westbesel/prediction/scripts/SCRATCH-1D_1.1/tmp/
RUN chmod ug+x /code/initialize.sh
# Expose ports
# 80 = Nginx
# 8000 = Gunicorn
# 3306 = MySQL
EXPOSE 80 8000 3306 2525

# Configure Nginx
RUN ln -s /code/nginx.conf /etc/nginx/sites-enabled/westbesel.conf
RUN rm /etc/nginx/sites-enabled/default

#Cron job
RUN apt-get install -y cron
RUN export EDITOR=nano
RUN touch /var/log/cron.log
COPY crontab /var/crontab.txt
RUN crontab /var/crontab.txt
RUN chmod 600 /etc/crontab

# Run Supervisor (i.e., start MySQL, Nginx, and Gunicorn)
COPY supervisord.conf /etc/supervisor/conf.d/supervisord.conf
CMD ["/usr/bin/supervisord"]
