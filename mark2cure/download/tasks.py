from __future__ import absolute_import

from .models import Download
from ..document.models import Document
from ..common import celery_app as app
from ..common.formatter import clean_df, apply_annotations

from celery import states
from celery.exceptions import SoftTimeLimitExceeded, Ignore
from django.utils import timezone

from django.core.files.storage import default_storage


@app.task(bind=True, ignore_result=False,
          max_retries=1, rate_limit='20/m', soft_time_limit=600,
          acks_late=True, track_started=True,
          expires=None)
def group_export(self, document_pks):
    try:
        save_location = 'downloads/data/bioc-{0}.xml'.format(timezone.now().strftime('%Y-%m-%d-%H-%M-%S'))

        writer = Document.objects.as_writer(documents=document_pks)
        org_er_df = Document.objects.entity_recognition_df(documents=document_pks, writer=writer)
        org_rel_df = Document.objects.relation_df(documents=document_pks)

        er_df = clean_df(org_er_df)
        export_writer = apply_annotations(writer, er_df=er_df, rel_df=org_rel_df)

        with default_storage.open(save_location, 'wb') as handle:
            handle.write(export_writer.__str__())

        d = Download.objects.create(
            # task_er=task_type == 'er',
            # task_rel=task_type == 'rel',
            task_er=True,
            task_rel=True,
            file=save_location
        )
        d.documents = document_pks
        d.save()

    except SoftTimeLimitExceeded:
        # (TODO) Retry?
        return

    if not self.request.called_directly:
        self.update_state(state=states.SUCCESS, meta=export_writer)
        raise Ignore()

