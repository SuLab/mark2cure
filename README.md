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
    . /var/www/virtualenvs/mark2cure-pro/bin/activate

    # Install all the python related dependencies
    easy_install -U distribute
    pip install -r service/requirements.txt

### Database Migrations


### Server Configuration

### Cron


### Control

`. /opt/mark2cure-venv/bin/activate`
`cd webapps/mark2cure/ && git pull`
`sudo supervisorctl restart mark2cure`

`python manage.py goldenmaster`

#### Notes

Print out a diagram of the database relationships: `python manage.py graph_models -a -o myapp_models.png`

Debug what anns are selected
YPet['730'].currentView.collection.parentDocument.get('annotations').each(function(m){ console.log(m.get('text')) })


from mark2cure.document.utils import *
user = User.objects.filter(username = 'demo').first()
doc = Document.objects.filter(pk = 375).first()
gold_matches(user,doc)



