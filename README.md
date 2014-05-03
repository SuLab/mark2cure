# Mark2Cure

[Mark2Cure](http://mark2cure.org/) is a web application concept recognition, relationship extraction and validation tool to generate annotations in fun and exploratory way to help scientific research!

Scientific communication is broken.

Progress in biomedical science is all about incrementally building off of the work of others. Currently, it takes far too long for one scientist's work to reach all the other scientists that would benefit. Mark2Cure will make scientific communication more efficient for everyone, leading to faster discoveries and cures for disease.

To be successful, we need your help. Mark2Cure works by directly involving crowds of people just like you.


## Setup

### Server Dependencies

* `sudo apt-get update`
* `sudo apt-get upgrade`
* `sudo apt-get install build-essential python python-dev python-pip python-virtualenv libmysqlclient-dev git-core nginx rabbitmq-server`

### Project Setup

* Make the project folder and download the repo
* `cd /var/www/ && git clone https://USER@bitbucket.org/sulab/mark2cure.git && cd mark2cure`

* Make the python virtual environment and activate it
* `virtualenv mark2cure-venv`
* `. /var/www/virtualenvs/mark2cure-venv/bin/activate`

* Install all the python related dependencies
* `easy_install -U distribute`
* `pip install -r requirements.txt`

### Database Migrations

* `sudo /opt/mark2cure-venv/bin/python manage.py schemamigration APP --auto CHANGE_MESSAGE`
* `sudo /opt/mark2cure-venv/bin/python manage.py migrate APP`

### Control

* `. /opt/mark2cure-venv/bin/activate`
* `cd webapps/mark2cure/ && git pull origin HEAD`
* `sudo supervisorctl restart mark2cure`

* `python /opt/python/current/app/manage.py celeryd -v 2 -E -l INFO`
* `python /opt/python/current/app/manage.py celerybeat`

### Utils

* Flow diagram of the database relationships
* `python manage.py graph_models -a -o myapp_models.png`

### Supervisor

/etc/supervisor/conf.d/mark2cure.conf

  [program:mark2cure]
  command = /bin/gunicorn_start
  user = deploy
  stdout_logfile = /home/ubuntu/webapps/mark2cure/logs/gunicorn_supervisor.log
  redirect_stderr = true


### NGINX

/etc/nginx/site-available/default

  server {
    listen 80;
    listen 443 ssl;

    server_name mark2cure.org www.mark2cure.org *.mark2cure.org;
    #root /var/www;

    access_log  /var/log/nginx/mark2cure.log;
    ssl_certificate /etc/ssl/mark2cure/mark2cure.crt;
    ssl_certificate_key /etc/ssl/mark2cure/mark2cure.key;

    location /static/admin {
      autoindex on;
      alias /opt/mark2cure-venv/lib/python2.7/site-packages/django/contrib/admin/static/admin;
    }

    location /static/grappelli {
      autoindex on;
      alias /opt/mark2cure-venv/local/lib/python2.7/site-packages/grappelli/static/grappelli/;
    }

    location /static {
      autoindex on;
      alias /home/ubuntu/webapps/mark2cure/static/;
    }

    location / {
          proxy_pass_header Server;
          proxy_set_header Host $http_host;
          proxy_redirect off;
          proxy_set_header X-Real-IP $remote_addr;
          proxy_set_header X-Scheme $scheme;
          proxy_connect_timeout 10;
          proxy_read_timeout 10;
          proxy_pass http://localhost:8080/;
    }

  }

### GUNICORN

  #!/bin/bash

  NAME="mark2cure_app"                              # Name of the application
  LOGFILE=/home/ubuntu/webapps/mark2cure/log/gunicorn.log
  LOGDIR=$(dirname $LOGFILE)
  DJANGO_SETTINGS_MODULE=mark2cure.settings
  DJANGO_WSGI_MODULE=mark2cure.wsgi

  DJANGODIR=/home/ubuntu/webapps/mark2cure                      # Django project directory
  VENVDIR=/opt/mark2cure-venv/bin/activate                       # Virtual envioronment directory

  USER=deploy                                       # the user to run as
  GROUP=deploy                                      # the group to run as
  NUM_WORKERS=3                                     # how many worker processes should Gunicorn spawn

  echo "Starting $NAME as `whoami`"

  # Activate the virtual environment
  source $VENVDIR

  cd $DJANGODIR
  test -d $LOGDIR || mkdir -p $LOGDIR

  exec /opt/mark2cure-venv/bin/gunicorn ${DJANGO_WSGI_MODULE}:application \
    --name $NAME \
    --workers $NUM_WORKERS \
    --user=$USER --group=$GROUP \
    --log-level=debug \
    --bind=0.0.0.0:8080
    --log-file=$LOGFILE 2>>$LOGFILE

#### Notes

* Debug what anns are selected
* `YPet['730'].currentView.collection.parentDocument.get('annotations').each(function(m){ console.log(m.get('text')) })`

Doc count: 3,454

from mark2cure.document.models import *
from django.contrib.auth.models import User

doc = Document.objects.filter(pk = 278).first()
user = User.objects.filter(username = 'worker').first()

doc.update_views(user, 'cr', True)


view = View.objects.filter(user__username = 'worker', section__document_id = 278).latest()
views = View.objects.filter(user__username = 'worker', section__document_id = 278).all()
for v in views:
  print v.id, v.created, v.completed, v.experiment


var $el = $('div.paragraph#556');
var svgContainer = d3.select($el[0]).append("svg").attr("width", $el.width()).attr("height", $el.height());

from django.contrib.auth.models import User
user = User.objects.get(username="admin")
user.userprofile.score()

