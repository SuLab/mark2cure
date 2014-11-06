from django.conf import settings
from django.db.models import signals
from django.core.mail import send_mail


import logging
logger = logging.getLogger(__name__)

