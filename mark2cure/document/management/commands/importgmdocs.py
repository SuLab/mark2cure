from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from mark2cure.document.models import Document, Section, View, Annotation
from mark2cure.document.utils import import_golden_documents, annotate_golden_documents, randomly_make_validation_documents

class Command(BaseCommand):
    args = '<gm_type>'
    help = 'Import GM documents'

    def handle(self, *args, **options):
        import_golden_documents()
        #randomly_make_validation_documents() #ALREADY RAN ONCE ON PROD
        annotate_golden_documents()
        self.stdout.write('Completed')

