# System Setup

This document describes a relatively detailed way to set up Mark2Cure so it could run locally on your computer.

The whole guide was tested on a Ubuntu machine. Mac users may find some difference when following this guide. Remember, whenever you have any trouble, Google is your best friend.

## Requirements

This guide assumes you know elementary knowledge on Linux command (including Git). If you are having trouble, feel free to visit [this website](https://www.google.com).

You will need a Linux system to set it up. If you are using a Mac, you may want to skip the following section. If you are using Windows OS, however, it is highly recommended that you use a virtual environment that runs Ubuntu. When setting up the local Mark2Cure, make sure you have a stable network connection.

**If you can talk to Max right now in person, you may want to ask him what `ENTREZ_EMAIL` is.** If he doesn't know what you are asking, tell him you are trying to fetch documents from PubMed.

### Install Ubuntu (skip if you are using Mac)

#### What you need
* VirtualBox
* Ubuntu .ISO image. You can download the desktop LTS version from Ubuntu website

#### Using VirtualBox
 The following guide assumes you are using VirtualBox, and unless otherwise specified, you can use the default settings
prompted. **Please take note of the bold part below.**

1. Download and install VirtualBox. Open VirtualBox.
2. Create a new virtual machine by clicking "New". Choose "Linux" and "Ubuntu" and give it a name.
3. For a smooth experience, you want to **allocate as much memory (RAM) for your virtual machine as possible**. Recommended memory size is at least 50% of the memory size of the machine you are working. For example, if you use a machine with 4G of RAM, you may want to allocate between 2500MB and 3000MB for this virtual machine.
4. Create a virtual hard disk (VHD). And **allocate at least 16GB of hard disk memory** for it. Failure to do so may result insufficient memory to install all the dependencies required to run Mark2Cure.
5. For the first time you run it, you need to select the original image file of Ubuntu. When prompted, choose the image file you just downloaded. Click "Install Ubuntu" and follow the instructions to customize.
6. Keep in mind that this is a virtual machine, so when it mentioned "Erase disk" or something similar, it won't do any hard to you real machine. Get yourself a cup of water while Ubuntu is finishing set-up, or keep reading to know what to do next.

### Other recommended program(s)
* MySQL Workbench.

  It is highly recommended that you download this for the sake of easier database processing you will need later. The guide will use MySQL Workbench when testing database.

# Deploy development environment

Before actually working on Mark2Cure, we need to make sure every dependency is installed. If you have your own preference or favorite commands, you can install them now.

Also, make a directory somewhere, and clone this GitHub repo. In this section, `root` or `./` refers to the directory where `.git` is located. For example, if your `.git` is located under `~/repos/mark2cure`, `./mark2cure/`
below will mean `~/repos/mark2cure/mark2cure/`

Some commands below assume that you have them installed already, so if your machine doesn't have these programs installed already (like if you have just installed Ubuntu), you will see many prompts that ask you to install them. Do  so when appropriate. For a better experience, a program called [f**k](https://github.com/nvbn/thefuck) can be very useful. There is a detailed instruction on how to install following the link.


## Dependencies

Before running the following commands, navigate to the root of the repo. **When prompted to put the password when installing MySQL, you need to remember that as it will be used later**
```
$ sudo apt-get update
$ sudo apt-get upgrade
$ sudo apt-get install build-essential python python-dev python-pip python-virtualenv libmysqlclient-dev git-core nginx supervisor rabbitmq-server graphviz libgraphviz-dev pkg-config libncurses5-dev npm ruby-dev
$ sudo pip install -r requirements.txt
$ sudo pip install nltk
$ sudo npm install gulp-cli -g
$ sudo npm install gulp -D
$ sudo npm install gulp-compass gulp-if gulp-livereload gulp-clean-css gulp-csso gulp-sass gulp-rename tiny-lr segfault-handler
$ sudo npm install -g bower
$ sudo bower install
$ sudo gem update --system
$ sudo gem install compass
$ sudo ln -s /usr/bin/nodejs /usr/bin/node
```
## Local setup

* Create a file at  `./mark2cure/local_settings.py` and put the following lines into the file. This file is what you want to modify if you want to customize anything locally. Note that previously you were asked to remember the password when installing MySQL, and now you will need to put it into `{YOUR PASSWORD HERE}` below. Also, if you get `ENTREZ_EMAIL` from Max, you need to put that here as well.
```
import os.path

SECRET_KEY = "thiscanactuallybesomerandomstring"
DEBUG = True
ACCOUNT_EMAIL_VERIFICATION = 'none'
ENTREZ_EMAIL = 'me@example.com' # Put that email address here

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': 'mark2cure',
        'USER': 'root',
        'PASSWORD': '{YOUR PASSWORD HERE}',
        'HOST': 'localhost',
        'PORT': ''
    }
}

PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
```

* We would like css and JavaScript files to be preprocessed. To do so, you should run this command **everytime the repo is updated**
  ```
  $ gulp
  ```

* Setup database

  * Open MySQL Workbench, and click something like "Local instance" (there should be only one if you just installed it).
  * In the Query tab, run `CREATE DATABASE mark2cure;`. You should see a new database is created under "SCHEMAS" under the bottom left panel.
  * Go back to console. Run `python manage.py migrate`. You will see that many data are being migrated.

* Import the basic training by running `python manage.py loaddata fixtures/tasks.json`. You will see a line at the end of the message like below. This imports the training data the user has to complete before they can actually do any tasks on Mark2Cure.

  ```
  Installed 8 object(s) from 1 fixture(s)
  ```

## Running server

Now we have finished setting up the database, it is time to run the server locally. However, make sure you have checked the tables under scheme `mark2cure` in the database. If you are seeing many tables, go ahead.

* Run the server by typing `python manage.py runserver`. If nothing goes wrong, you should see something below.

  However, keep in mind that once you run this command, you may want to open a new terminal/console. Also, every time you modify something, make sure you stop and re-run the server.
  ```
  System check identified no issues (0 silenced).
  Django version 1.9.6, using settings 'mark2cure.settings
  Starting developement server at http://127.0.0.1:8000/
  Quit the server with CONTROL-C
  ```


* Visit `127.0.0.1:8000` in your broswer. You will see the front page of Mark2Cure. From here, you are running a local Mark2Cure. The experience should be exactly the same as you would on `mark2cure.org`

## Create a new account

Sign up a new account as you would on `mark2cure.org` by going over those training steps. Remember the username and password, as you will use this account most of time later. For the email, it does not actually matter what it is as long as it looks like an email address (since we have disabled it by `ACCOUNT_EMAIL_VERIFICATION = 'none'` above)

### Doing training tasks

You will need to do some training before you actually see the dashboard (as you would on the real Mark2Cure). To create groups and see those quests, you also need to finish all those tasks listed on the left.

### ... or not

If you really do not want to do them, here is an alternative that you could do to skip them by modifying the database. **However, these scripts are created by the author of this guide based on his sole observation of how Mark2Cure works**. The author has included his analysis and you can modify the SQL scripts accordingly before applying them.

Mark2Cure stores the progress of each user's tasks in the table called `task_level`. There are different levels for different tasks, so you can append rows here to skip them, for example, by running the following scripts. If you are not very sure what you are doing, also read the note below before applying it.

  ```
  INSERT INTO `mark2cure`.`task_level` (`id`, `task_type`, `level`, `created`, `user_id`) VALUES ('1000000', 'ner', '3', '1999-09-09 09:09:09.999999', '1');
  INSERT INTO `mark2cure`.`task_level` (`id`, `task_type`, `level`, `created`, `user_id`) VALUES ('1000001', 're', '1', '1999-09-09 09:09:09.999999', '1');
  INSERT INTO `mark2cure`.`task_level` (`id`, `task_type`, `level`, `created`, `user_id`) VALUES ('1000002', 're', '2', '1999-09-09 09:09:09.999999', '1');
  INSERT INTO `mark2cure`.`task_level` (`id`, `task_type`, `level`, `created`, `user_id`) VALUES ('1000003', 're', '3', '1999-09-09 09:09:09.999999', '1');
  INSERT INTO `mark2cure`.`task_level` (`id`, `task_type`, `level`, `created`, `user_id`) VALUES ('1000004', 'ner', '4', '1999-09-09 09:09:09.999999', '1');
  INSERT INTO `mark2cure`.`task_level` (`id`, `task_type`, `level`, `created`, `user_id`) VALUES ('1000005', 'ner', '5', '1999-09-09 09:09:09.999999', '1');
  INSERT INTO `mark2cure`.`task_level` (`id`, `task_type`, `level`, `created`, `user_id`) VALUES ('1000006', 'ner', '6', '1999-09-09 09:09:09.999999', '1');
  INSERT INTO `mark2cure`.`task_level` (`id`, `task_type`, `level`, `created`, `user_id`) VALUES ('1000007', 'ner', '7', '1999-09-09 09:09:09.999999', '1');
  ```

Two things that you need to know before executing:

1. In the tuple after `VALUES` in each line, the first element (e.g. `'1000000'`, `'1000001'`, etc. has to be unique.
2. You need to change `user_id` if you are not applying it for the first user created (in the snippet above, `1` was used). If you need to know the `user_id` of the account you just created, view the content in the table called `account_emailaddress`.

## Create groups of Quests

If you do not know what a Quest is, Quest is a set of Documents you are going to annotate to identify possible concepts (like gene, drug or treatment). The data will be later used in relation tasks.

### Fetch documents from PubMed

We first need to fetch the documents from PubMed, and then we will use these documents to create quests.

Documents in PubMed are identified by a unique ID, which is a string of digits. If you want to import specific documents from PubMed, you may want to write them down.

Now, run `python manage.py shell`, which will open a console and place the cursor after `In [1]: `.

#### First time setup

If this is NOT the first time trying to fetch documents from PubMed, you can skip this part.

You need to download `nltk` library before fetching any documents. To do so, in the shell opened, type in

```
import nltk
nltk.download()
```

This will open up an NLTK Downloader, which prompts you five options. Enter `d` to download, and put `all` when asked which package to download. The downloading may take some time.

#### Batch fetch or Fetch a specific document

Type the following scripts into the shell, or you can modify it to fetch certain documents. Please note that there might be a minimum number of documents required to create a quest

```
import mark2cure.document.tasks

# Importing some random documents
for x in range (27834000, 27834100):
  mark2cure.document.tasks.get_pubmed_document(x);

# Or, you can import a specific document given an ID:
mark2cure.document.tasks.get_pubmed_document(27834101);
```

### Create groups

After the documents are fetched, type the following into the shell. Note that you can edit the group property (name and description) later through the database.

```
from mark2cure.common.models import Group, Document

G = Group()
G.save()
G.assign(Document.objects.all()) # Here, you can modify it to only put specific documents
G.save()
G.enabled = True
G.save()
```

To know more about `assign` function used above, read the source code [here](https://github.com/SuLab/mark2cure/blob/13b6170328b4fff0ca720450eac5e8b14ecc4ad7/mark2cure/common/models.py#L126) (If the link is broken, refer to `assign` function in `./common/models.py`).
.

# Credit

The original document was written by [Runjie Guan](https://github.com/AnoXDD), who was an intern working for this project between Sep 2016 to Mar 2017. I have gone over the same procedure described above and can make sure it is accurate as of Mar 2017. This guide is not comprehensive since I did not have the chance to work on every feature of Mark2Cure, so please contact Max if you need help. For example, this guide doesn't tell people how to create relationship tasks, so if you do, please add it into this guide to complete it.

Since the project is developing quite quickly, if you see any problems in this guide, you will probably want to figure it out by yourself, which is what I did before I know how to write this.




========


## Setup

### Server Dependencies

* `sudo apt-get update`
* `sudo apt-get upgrade`
* `sudo apt-get install build-essential python python-dev python-pip python-virtualenv libmysqlclient-dev git-core nginx supervisor rabbitmq-server graphviz libgraphviz-dev pkg-config libncurses5-dev`

### Project Setup

* Make the python virtual environment and activate it
* `virtualenv mark2cure-venv`
* `. /var/www/virtualenvs/mark2cure-venv/bin/activate`

* Make the project folder and download the repo
* `sudo adduser deploy`
* `sudo /home/deploy/webapps`
* `cd /home/deploy/webapps/`
* `git clone https://USER@bitbucket.org/sulab/mark2cure.git && cd mark2cure`

* Install all the python related dependencies
* `sudo /opt/mark2cure-venv/bin/pip install -r requirements.txt`


### Database Migrations

* `python manage.py schemamigration APP --auto CHANGE_MESSAGE`
* `python manage.py migrate APP`

### Control

* `. /opt/mark2cure-venv/bin/activate`
* `cd webapps/mark2cure/ && git pull origin HEAD`
* `sudo supervisorctl restart mark2cure`

* `sudo chmod a+x /bin/gunicorn_start`

### Utils

* Flow diagram of the database relationships
* `python manage.py graph_models -a -o myapp_models.png`


#### Notes

* Networkx talk: http://www.slideshare.net/arnicas/a-quick-and-dirty-intro-to-networkx-and-d3



# **Mark2Cure** detailed installation instructions for developers

These are instructions, tips, and tricks for how to install **mark2cure** for
developers with more detail than above, assuming you might be unfamiliar
with much of this technology. Please feel free to let us know if anything is
unclear.

Below is the terminal history for installation of **mark2cure** in a
“clean” brand new Macbook Pro. You will likely have a computer that has already
been used for development, so you may run into different issues, but hopefully
none.

First thing!
You will need pip (you probably already have this). pip is important because it
is compatible with virtual environments which we will use to compartmentalize
different projects that require different versions of software. Mark2cure has a
lot of dependencies.

Install [Brew] using the long command below. It’s a package manager for Apple
comps. “installs the stuff you need that Apple didn’t.”

```sh
$ ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Sudo install a few things. You will need to have virtual environments
installed. For more information, see [Python Virtual Environments].

```sh
$ sudo pip install virtualenv
```

You will also need the virtual environment wrapper: [Virtual Environment Wrapper]

```sh
$ sudo pip install virtualenvwrapper
```

Open the /.bash_profile.

```sh
$ nano ~/.bash_profile
```

Add the following *two* lines to your bash profile:

export WORKON_HOME=$HOME/.virtualenvs

source /usr/local/bin/virtualenvwrapper.sh

Then, "source" the profile.

```sh
$ source ~/.bash_profile
```

Make a new directory where you will house your development
version of **mark2cure** (and other development projects you may have).
You could call this something like "repos." Go into this directory.

```sh
$ mkdir repos
$ cd repos/
```

Clone the repository so that you have an exact copy of **mark2cure**. This would
be done in your “repos” folder. No need to make a folder called “mark2cure”
because “git clone” will do this for you. If you want to develop, then you would
make a branch later (see git documentation).

```sh
$ git clone https://your_user_name@bitbucket.org/sulab/mark2cure.git
```

cd into the newly created folder called **mark2cure**, and make a
virtual environment called “**mark2cure**."

```sh
$ mkvirtualenv mark2cure
```

You want to *activate* the virtual environment.

```sh
$ workon mark2cure
```

There are two folders called **mark2cure**. One is a higher level folder, but
you should "cd" into the lower level to install the requirements file. The
requirements file is just a list of software programs that **mark2cure** (or
any program with a requirements.txt file) is dependent on.

```sh
$ cd mark2cure/
$ pip install -r requirements.txt
```

Databases are needed in **mark2cure** to save and recall all of the useful
annotations provided by the Mark2Curators.

Make sure you download MySQL [Sequel Pro] and obtain a copy of the MySQL
"framework" **mark2cure** database. See Max or Jennifer for this.

```sh
$ brew install mysql
$ ln -sfv /usr/local/opt/mysql/*.plist ~/Library/LaunchAgents
$ launchctl load ~/Library/LaunchAgents/homebrew.mxcl.mysql.plist
```

If you have issues with graphviz and pygraphviz, you can remove them from the
requirements.txt file (*for now*). Other dependencies probably should not be
tampered with.

```sh
$ cd repos/
$ cd mark2cure/
$ pip install -r requirements.txt
```

Make new directory called “env_vars”

```sh
$ mk_dir env_vars
$ mkdir env_vars
$ cd env_vars/
$ touch development.sh
$ nano development.sh
```

Inside the development.sh file will live some special information.
Please see Max or Jennifer for the file contents.

```sh
$ cd mark2cure/
$ source env_vars/development.sh
$ echo $MARK2CURE_DATABASE_URL
```

You will also need the NLTK database [Natural Language Tool Kit]. Download NLTK
using the following commands and follow download instructions in the popup
window; download *everything*.

```sh
$ sudo pip install -U nltk
$ python
```
Running the Python shell will allow you test that the download was successful:

```python
import nltk
```

Run the unit tests:

```sh
$ python manage.py test
```

At this point, all tests should be passing. After all the tests
work correctly, then you can run the server.

```sh
$ python manage.py runserver_plus
```

You should now have a *local* development version of Mark2Cure!

Install “tree” so that you can see a data tree of your current working
directory to see the Mark2Cure application layout.

```sh
$ brew install tree
```

# **mark2cure** development details

When doing manual testing of the **mark2cure** application (i.e., using the
website, and not using the command `python manage.py test`), you will go to the
development server site([Django development server]), and you will realize that
you cannot get past the "log in" page. To fix this, you can run these commands
to *force* the development server to let you get to the "test_user" account.

From terminal, use `python manage.py shell_plus` and run the following
Python commands:

```python
User.objects.create_user('test_user', password='password') 
```

This does not solve the issue of a "simulated database with real content for
testing"... more to come.

## Misc. Information about the project (to be expanded):

[Ypet] is a Javascript library built in Marionette.js to rapidly annotate
paragraphs of text on websites. YPet was developed to rapidly annotate
bio-medical literature for **mark2cure** at The Su Lab by Max Nanis.

## Troubleshooting:

You will need to create an account to [Sentry] to log in and monitor error codes
and issues with both development and production versions of **mark2cure**.

When running `$ python manage.py runserver_plus"` did you get this error?:

```sh
Secret value 'SENTRY_PUBLIC_KEY' is not set
```

This indicates that you must do the following:

```sh
$ echo $MARK2CURE_SENTRY_PUBLIC_KEY
```
Oh no, nothing is there, so you must "source" the development file:

```sh
$ source env_vars/development.sh
```

```sh
$ echo $MARK2CURE_SENTRY_PUBLIC_KEY
```
There, now you should see the key output and can run
`$ python manage.py runserver_plus"` successfully for testing on the local
server.

[Virtual Environment Wrapper]:http://virtualenvwrapper.readthedocs.org/en/latest/install.html
[Django Web Framework]:https://www.djangoproject.com/start/overview/
[Python Virtual Environments]:http://docs.python-guide.org/en/latest/dev/virtualenvs/
[Ypet]:https://github.com/SuLab/YPet
[Dillinger]:http://dillinger.io/
[Sentry]:http://sentry.sulab.org/tsri/development/
[Brew]:http://brew.sh/
[Sequel Pro]:http://www.sequelpro.com/
[Natural Language Tool Kit]:http://www.nltk.org/
[Mark2Cure.org]:https://mark2cure.org/
[Django development server]:http://127.0.0.1:8000/

## Migration Commands

------ training requirement migrations -----------

for lvl in Level.objects.filter(task_type='re'):
  r = Requirement.objects.filter(task_type='re', order=lvl.level).first()
  if(r):
    lvl.requirement = r
    lvl.save()


for u in User.objects.all():
  for lvl in Level.objects.filter(task_type='ner', user_id=u.pk).order_by('-level'):

      r = Requirement.objects.filter(task_type='ner', order=lvl.level+1).first()
      if(r):
        lvl.requirement = r
        lvl.save()

