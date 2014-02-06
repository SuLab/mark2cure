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
`cd webapps/mark2cure/ && git pull origin HEAD`
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

In [1]: from mark2cure.document.utils import *

In [2]: user = User.objects.filter(username = 'demo').first()

In [3]: doc = Document.objects.filter(pk = 375).first()


In [4]: gold_matches(user,doc)
user_annotations: []
 ~ ~ ~ ~ ~
gold_annotations: [u’retinoblastoma', u'retinoblastoma']
 ~ ~ ~ ~ ~
true_positives: []
 ~ ~ ~ ~ ~
Out[4]: 0

In [5]: gold_matches(user,doc)
user_annotations: [u’retinoblastoma gene']
 ~ ~ ~ ~ ~
gold_annotations: [u’retinoblastoma', u'retinoblastoma']
 ~ ~ ~ ~ ~
true_positives: []
 ~ ~ ~ ~ ~
Out[5]: 0

In [6]: gold_matches(user,doc)
user_annotations: [u’retinoblastoma gene', u'retinoblastoma']
 ~ ~ ~ ~ ~
user_annotations: [u’retinoblastoma', u'retinoblastoma']
 ~ ~ ~ ~ ~
true_positives: [u’retinoblastoma', u'retinoblastoma']
 ~ ~ ~ ~ ~
Out[6]: 2



# Relationship Builder

## Open Questions
  * Does the user make the combinations or do they confirm combinations

## Potential thought mechanics
  * http://en.wikipedia.org/wiki/Concentration_(game)
  * http://connectmania.com/

## Tools
  * http://www.mindmeister.com/

    I like the link in the node area, click on it to take you to an external URL, but could link to a different concept_url

    Relationship dragging sucks. http://www.puff.me.uk/ss/JECWjv6O3OTz3.png

    $6 / $10 / $15 a month

  * http://dev.slatebox.com/

    Waste of time.

## Sketchs

  * http://jqueryui.com/droppable/#revert
  * http://www.puff.me.uk/ss/6NWJZywgEeLv0.png

  a,b,c,d are divs that are all draggable — when you drag it becomes lighter and if you drop it on another concept div, it creates a list of relationships else where on the page (below the line). makes the assumption that a concept can’t be related to itself but that’s okay (no b<==>b groups)
the “relationships” would have some controls like delete, etc to change after they’ve been made


how others do drill down selection and help definitions
