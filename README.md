## Setup

### Local

. ENV/bin/activate

### Server

#### Dependencies

sudo apt-get install nodejs

#### Setup

The application runs with uWSGI. To deploy:

git clone git@bitbucket.org:sulab/mark2cure.git /var/www/
Add a MySQL DB connection details
pip install -r service/requirements.txt
/usr/local/bin/uwsgi --ini /var/www/mark2cure/deploy/uwsgi.ini --daemonize /var/www/mark2cure/deploy/uwsgi.log

### Cron

*/30 * * * * /var/www/abstract-annotations/ENV/bin/python /var/www/abstract-annotations/service/manage.py heatmap
0 */2 * * * /var/www/abstract-annotations/ENV/bin/python /var/www/abstract-annotations/service/manage.py annotate
