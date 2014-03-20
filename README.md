# Mark2Cure

[Mark2Cure](http://mark2cure.org/) is a web application concept recognition, relationship extraction and validation tool to generate annotations in fun and exploratory way to help scientific research!

Scientific communication is broken.

Progress in biomedical science is all about incrementally building off of the work of others. Currently, it takes far too long for one scientist's work to reach all the other scientists that would benefit. Mark2Cure will make scientific communication more efficient for everyone, leading to faster discoveries and cures for disease.

To be successful, we need your help. Mark2Cure works by directly involving crowds of people just like you.


## Setup

### Server Dependencies

* `sudo apt-get update`
* `sudo apt-get upgrate`
* `sudo apt-get install build-essential python python-dev python-pip python-virtualenv libmysqlclient-dev git-core nginx`

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

### Utils

* Flow diagram of the database relationships
* `python manage.py graph_models -a -o myapp_models.png`

#### Notes

* Debug what anns are selected
* `YPet['730'].currentView.collection.parentDocument.get('annotations').each(function(m){ console.log(m.get('text')) })`


