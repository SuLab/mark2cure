from django.conf import settings
from django.db.models import Count
from django.core.mail import send_mail

from mark2cure.document.models import Activity

from random import shuffle
import logging
logger = logging.getLogger(__name__)


