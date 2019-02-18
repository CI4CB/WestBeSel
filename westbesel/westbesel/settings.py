"""
Django settings for django_docker project.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.6/ref/settings/
"""

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
BASE_DIR = os.path.dirname(os.path.dirname(__file__))


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.6/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'w5-iqfj^5cl6tj0_8*xmgnj5r^@&p@n+)=6g1-!^i0&sohyicg'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = False

ALLOWED_HOSTS = ['*']

# Application definition

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'prediction',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'westbesel.urls'

WSGI_APPLICATION = 'westbesel.wsgi.application'


# Database
# https://docs.djangoproject.com/en/1.6/ref/settings/#databases

#DATABASES = {
#    'default': {
#        'ENGINE': 'django.db.backends.mysql', # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
#        'NAME': 'devdb',                      # Or path to database file if using sqlite3.
#        # The following settings are not used with sqlite3:
#        'USER': 'devuser',
#        'PASSWORD': 'devpass', # Entered via fab command; leave blank if using SQLite
#        'HOST': '',                      # Empty for localhost through domain sockets or '127.0.0.1' for localhost through TCP.
#        'PORT': '',                      # Set to empty string for default.
#    }
#}


DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
        #'NAME': '/code/westbesel/db/db.sqlite3',
    }
}

#DATABASES = {
#  'default': {
#    'ENGINE': 'django.db.backends.mysql', 
#    'NAME': 'djangodocker_db',
#    'USER': "root",
#    'PASSWORD': 'root',
#    'HOST': '127.0.0.1', # Or an IP Address that your DB is hosted on
#    'PORT': '3306',
#  }
#}



# Internationalization
# https://docs.djangoproject.com/en/1.6/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'Europe/Paris'

USE_I18N = True

USE_L10N = True

USE_TZ = False


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.6/howto/static-files/

STATIC_URL = '/static/'

STATICFILES_DIRS = (
    os.path.join(BASE_DIR, "static"),
)

STATIC_ROOT = '/static/';

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            os.path.join(BASE_DIR, 'templates'),
        ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                # Insert your TEMPLATE_CONTEXT_PROCESSORS here or use this
                # list if you haven't customized them:
                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.i18n',
                'django.template.context_processors.media',
                'django.template.context_processors.static',
                'django.template.context_processors.tz',
                'django.contrib.messages.context_processors.messages',
            ],
            'debug': True,
            'libraries':{
            'mytags': 'prediction.templatetags.mytags',

            },
        },
    },
]

# Import production settings if the environment variable DJANGO_PRODUCTION is true
# (It's set to True by default in the Dockerfile, but you can override it with `docker run` for development)
#if os.environ['DJANGO_PRODUCTION'] == 'true':
#    from settings_production import *
