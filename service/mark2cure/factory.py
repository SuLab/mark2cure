# -*- coding: utf-8 -*-
"""
    mark2cure.factory
    ~~~~~~~~~~~~~~~~

    mark2cure factory module
"""

import os

from celery import Celery
from flask import Flask

from .core import db, mail, login_manager
from .models import User
from .middleware import HTTPMethodOverrideMiddleware

def create_app(package_name, settings_override=None):
    """Returns a :class:`Flask` application instance configured with common
    functionality for the Mark2Cure platform.

    :param package_name: application package name
    :param settings_override: a dictionary of settings to override
    """
    app = Flask(package_name, instance_relative_config=True,
        static_url_path = '',
        static_folder = '../../../web-app')

    app.config.from_object('mark2cure.settings')
    app.config.from_object(settings_override)

    db.init_app(app)
    mail.init_app(app)
    login_manager.init_app(app)

    app.wsgi_app = HTTPMethodOverrideMiddleware(app.wsgi_app)

    return app


def create_celery_app(app=None):
    app = app or create_app('mark2cure', os.path.dirname(__file__))
    celery = Celery(__name__, broker=app.config['CELERY_BROKER_URL'])
    celery.conf.update(app.config)
    TaskBase = celery.Task

    class ContextTask(TaskBase):
        abstract = True

        def __call__(self, *args, **kwargs):
            with app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)

    celery.Task = ContextTask
    return celery

