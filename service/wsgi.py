# -*- coding: utf-8 -*-
"""
    wsgi
    ~~~~

    mark2cure wsgi module
"""

from werkzeug.serving import run_simple
from werkzeug.wsgi import DispatcherMiddleware

from mark2cure import api
from mark2cure.core import login_manager

# Create user loader function
@login_manager.user_loader
def load_user(user_id):
    print "/ / / / / / / / / loading user / / / / / / / / / /"
    return db.session.query(User).get(user_id)


application = DispatcherMiddleware(api.create_app())

if __name__ == "__main__":
    run_simple('0.0.0.0', 5000, application, use_reloader=True, use_debugger=True)
