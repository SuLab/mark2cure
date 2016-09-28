from __future__ import unicode_literals

from django.apps import AppConfig


class DocumentConfig(AppConfig):
    name = 'mark2cure.document'
    verbose_name = 'Pubmed Document'

    def ready(self):
        import mark2cure.document.signals # noqa
