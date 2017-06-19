from __future__ import absolute_import

from .models import Download
from ..document.models import Document
from ..common.formatter import clean_df, apply_annotations

# from celery import states
# from celery.exceptions import SoftTimeLimitExceeded
# from ..common import celery_app as app

from django.core.files.storage import default_storage

from django.utils import timezone
import random


# @app.task(bind=True, ignore_result=True,
#           max_retries=1, rate_limit='20/m', soft_time_limit=600,
#           acks_late=True, track_started=True,
#           expires=None)
def group_export(self, document_pks, export_type=0):
    pass
#     save_location = 'downloads/data/bioc-{0}.xml'.format(timezone.now().strftime('%Y-%m-%d-%H-%M-%S'))
#     er_df = None
#     rel_df = None
#
#     try:
#         # Base BioC writer for all exports
#         writer = Document.objects.as_writer(documents=document_pks)
#
#         if export_type == 0 or export_type == 1:
#             er_df = clean_df(Document.objects.entity_recognition_df(documents=document_pks, writer=writer))
#
#         if export_type == 0 or export_type == 2:
#             rel_df = Document.objects.relation_df(documents=document_pks)
#
#         export_writer = apply_annotations(writer, er_df=er_df, rel_df=rel_df)
#
#         with default_storage.open(save_location, 'wb') as handle:
#             handle.write(export_writer.__str__())
#
#         d = Download.objects.create(
#             task_er=True if export_type == 0 or export_type == 1 else False,
#             task_rel=True if export_type == 0 or export_type == 2 else False,
#             file=save_location
#         )
#         d.documents = document_pks
#         d.save()
#
#     except SoftTimeLimitExceeded:
#         return False
#     except:
#         self.retry(countdown=int(random.uniform(2, 4) ** self.request.retries))
#
#     if not self.request.called_directly:
#         self.update_state(state=states.SUCCESS, meta=export_writer)
#         return True
#
