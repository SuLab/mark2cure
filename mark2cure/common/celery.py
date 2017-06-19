# from __future__ import absolute_import
#
# from django.conf import settings
#
# from celery import Celery, states
# from celery.exceptions import Ignore
#
# import logging
# import os
#
# logger = logging.getLogger(__name__)
#
# # set the default Django settings module for the 'celery' program.
# os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'mark2cure.settings')
#
# app = Celery('mark2cure')
#
# # Using a string here means the worker will not have to
# # pickle the object when using Windows.
# app.config_from_object('django.conf:settings')
# app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)
#
#
# @app.task(bind=True)
# def debug_task(self, msg):
#     logger.debug('Debug Task', exc_info=True, extra={'self': self})
#
#     if not self.request.called_directly:
#         self.update_state(state=states.SUCCESS, meta=msg)
#         raise Ignore()
#
