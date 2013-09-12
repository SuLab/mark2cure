# Mark2Cure

An online text annotator to generate user annotations in fun and exploratory way

## Setup (Ubuntu Server 12.04)

### Server Dependencies

    sudo apt-get update
    sudo apt-get upgrate

    sudo apt-get install build-essential python python-dev python-pip python-virtualenv libmysqlclient-dev
    sudo apt-get install git-core nginx nodejs

### Project Setup

    # make the project folder and download the repo
    cd /var/www/ && git clone https://USER@bitbucket.org/sulab/mark2cure.git && cd mark2cure

    # Make the python virtual environment and activate it
    virtualenv ENV
    . ENV/bin/activate

    # Install all the python related dependencies 
    easy_install -U distribute
    pip install -r service/requirements.txt

    # Start the server to start accepting connections
    uwsgi --ini /var/www/mark2cure/deploy/uwsgi.ini --daemonize /var/www/mark2cure/deploy/uwsgi.log

### Server Configuration

    # /etc/nginx/sites-available/default
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

    */30 * * * * /var/www/mark2cure/ENV/bin/python /var/www/mark2cure/service/manage.py heatmap
    0 */2 * * * /var/www/mark2cure/ENV/bin/python /var/www/mark2cure/service/manage.py annotate

### Control

    # View the server log
    tail -f /var/www/mark2cure/deploy/uwsgi.log
    # Restart gracefully
    kill -HUP $(cat /var/www/mark2cure/deploy/uwsgi_master.pid)
    # Kill
    kill -INT $(cat /var/www/mark2cure/deploy/uwsgi_master.pid)

### Update

    #!/bin/bash
    LOGPATH=/home/ubuntu/gitpull.log

    echo "Updating repo">>$LOGPATH
    cd /var/www/mark2cure && git pull origin master

    echo "Compiling Application">>$LOGPATH
    cd /var/www/mark2cure/web-app && node r.js -o app.build.js

    echo "Restarting uWSGI">>$LOGPATH
    kill -HUP $(cat /var/www/mark2cure/deploy/uwsgi_master.pid)

    date "+%F %T">>$LOGPATH


#### Notes

db = SQLAlchemy(app)
admin = db.session.query(User).get(1)
doc = db.session.query(Document).get(28)

quest = Quest('super interesting list', admin)
db.session.add(quest)
db.session.commit()

qr = QuestRelation(quest, db.session.query(Document).get(29))
db.session.add(qr)
db.session.commit()
