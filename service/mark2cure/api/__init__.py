# -*- coding: utf-8 -*-
"""
    mark2cure.api
    ~~~~~~~~~~~~~

    mark2cure api application package
"""

from flask.ext.restful import Api

from ..core import db
from ..core import Mark2CureError, Mark2CureFormError
from ..helpers import JSONEncoder
from .. import factory

from ..routes import *

# Import the api endpoints
from .annotations import Annotations
from .documents import Documents
from .messages import Messages
from .network import Network
from .users import Users
from .email import Email
from .gold import Gold
from .quests import Quests

def create_app(settings_override=None):
    """Returns the Mark2Cure API application instance"""

    app = factory.create_app(__name__, settings_override)

    # Set the default JSON encoder
    # app.json_encoder = JSONEncoder

    # Register custom error handlers
    # app.errorhandler(Mark2CureError)(on_overholt_error)
    # app.errorhandler(Mark2CureFormError)(on_overholt_form_error)
    # app.errorhandler(404)(on_404)


    api = Api(app)

    api.add_resource(Email,       '/api/v1/email')
    api.add_resource(Users,       '/api/v1/user')
    api.add_resource(Messages,    '/api/v1/messages')
    api.add_resource(Annotations, '/api/v1/annotations',  '/api/v1/annotations/<int:ann_id>')
    api.add_resource(Documents,   '/api/v1/documents',    '/api/v1/documents/<int:doc_id>')
    api.add_resource(Network,     '/api/v1/network')
    api.add_resource(Gold,        '/api/v1/gm',           '/api/v1/gm/<int:doc_id>')
    api.add_resource(Quests,      '/api/v1/quest',        '/api/v1/quest/<int:quest_id>',
                                                          '/api/v1/quest/<string:quest_name>')

    return app


