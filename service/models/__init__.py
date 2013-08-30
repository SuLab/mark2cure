from flask.ext.sqlalchemy import SQLAlchemy

db = SQLAlchemy()

from .annotation import Annotation
from .concept import Concept
from .document import Document
from .message import Message
from .user import User
from .view import View
