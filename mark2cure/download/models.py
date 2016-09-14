from django.db import models


class Download(models.Model):
    documents = models.ManyToManyField('document.Document')

    task_er = models.BooleanField(default=False)
    task_rel = models.BooleanField(default=False)

    file = models.FileField(null=True, blank=True)
    create_time = models.DateTimeField(auto_now_add=True)

    download_count = models.IntegerField(default=0)

    class Meta:
        app_label = 'download'

