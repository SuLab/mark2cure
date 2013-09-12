from flask import Flask
from flask.ext.restful import Api

from app import app
from api import *
from models import *

import IPython, settings

app.testing = True
test_client = app.test_client()

welcome_message = """Welcome to your Flask CLI environment. 
The following variables are available to use:

app           -> Your Flask app instance.
test_client   -> Your Flask app.test_client().
"""

IPython.embed(header=welcome_message)


#         user = User('username',
#                     'email',
#                     1,
#                     'asdf')
#         db.session.add(user)
#         db.session.commit()
# 
# db = SQLAlchemy()
# def create_app():
#     app = Flask(__name__)
#     db.init_app(app)
#     with app.test_request_context():
#         db.create_all()
#     return app
