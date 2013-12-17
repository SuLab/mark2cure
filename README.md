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

. /var/www/virtualenvs/mark2cure-prod/bin/activate
sudo python /var/www/mark2cure/manage.py run_gunicorn -w 4 -k gevent

#### Notes

Print out a diagram of the database relationships: `python manage.py graph_models -a -o myapp_models.png`
