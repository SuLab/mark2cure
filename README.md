## Setup

### Server

#### Dependencies

sudo apt-get update
sudo apt-get upgrate

sudo apt-get install build-essential python python-dev python-pip python-virtualenv libmysqlclient-dev
sudo apt-get install git-core uwsgi uwsgi-plugin-python nginx nodejs

pip install uwsgi

#### Setup

The application runs with uWSGI. To deploy:

  cd /var/www/
  git clone git@bitbucket.org:sulab/mark2cure.git && cd mark2cure
  mark2cure
  virtualenv ENV
  . ENV/bin/activate
  easy_install -U distribute
  pip install -r service/requirements.txt

  /usr/bin/uwsgi --ini /var/www/mark2cure/deploy/uwsgi.ini --daemonize /var/www/mark2cure/deploy/uwsgi.log

#### Config

  server {
    listen 80;

    root /var/www/mark2cure/web-app;
    index index.html;

    server_name mark2cure.org;

    location /api {
      include uwsgi_params;
      uwsgi_pass 127.0.0.1:3031;
    }
  }

### Cron

*/30 * * * * /var/www/abstract-annotations/ENV/bin/python /var/www/abstract-annotations/service/manage.py heatmap
0 */2 * * * /var/www/abstract-annotations/ENV/bin/python /var/www/abstract-annotations/service/manage.py annotate
