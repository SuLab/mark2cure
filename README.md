## Setup

### Server

#### Dependencies

sudo apt-get update
sudo apt-get upgrate
sudo apt-get install build-essential python python-dev python-pip git-core uwsgi nginx nodejs
pip install uwsgi virtualenv

#### Setup

The application runs with uWSGI. To deploy:

git clone git@bitbucket.org:sulab/mark2cure.git /var/www/
Add a MySQL DB connection details
pip install -r service/requirements.txt
/usr/bin/uwsgi --ini /var/www/mark2cure/deploy/uwsgi.ini --daemonize /var/www/mark2cure/deploy/uwsgi.log

### Cron

*/30 * * * * /var/www/abstract-annotations/ENV/bin/python /var/www/abstract-annotations/service/manage.py heatmap
0 */2 * * * /var/www/abstract-annotations/ENV/bin/python /var/www/abstract-annotations/service/manage.py annotate
