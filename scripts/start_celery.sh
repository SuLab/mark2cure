#!/bin/sh

# Start the venv
source /opt/mark2cure-venv/bin/activate
source /opt/mark2cure-venv/bin/activate

# Make sure rabbit-mq is running and start celery worker and beat
sudo /etc/init.d/rabbitmq-server start
python /home/ubuntu/webapps/mark2cure/manage.py celeryd_detach worker --pidfile=/home/ubuntu/webapps/mark2cure/logs/celery_worker.pid
python /home/ubuntu/webapps/mark2cure/manage.py celeryd_detach beat --pidfile=/home/ubuntu/webapps/mark2cure/logs/celery_beat.pid

