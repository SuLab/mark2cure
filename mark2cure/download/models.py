from django.db import models
import os

app_path = 'download/'


def _download_file_name(instance, filename, folder=app_path + 'data/bioc'):
    name = os.urandom(40).encode('hex') + os.path.splitext(filename)[1]
    return '/'.join([folder, name])


class Download(models.Model):
    documents = models.ManyToManyField('document.Document')

    task_er = models.BooleanField(default=False)
    task_rel = models.BooleanField(default=False)

    file = models.FileField(upload_to=_download_file_name, null=True, blank=True)

    create_time = models.DateTimeField(auto_now_add=True)

    download_count = models.IntegerField(default=0)

    class Meta:
        app_label = 'download'

