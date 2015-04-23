from django.core.management.base import BaseCommand
from optparse import make_option
from django.contrib.auth.models import User
from django.conf import settings

from mark2cure.document.models import Document, Pubtator, Section, View, Annotation

from mark2cure.common.bioc import BioCReader, BioCDocument, BioCPassage

from collections import Counter
import requests
import random
import csv


class Command(BaseCommand):
    help = 'Utils for populating the system with Documents & Annotations'

    option_list = BaseCommand.option_list + (
        make_option('--keys',
            action='store_true',
            dest='keys',
            default=False,
            help='Inspect the bioc infon keys'),
    )

    def handle(self, *args, **options):
        types_arr = []
        errors = 0

        if options['keys']:
            for pubtator in Pubtator.objects.filter(content__isnull=False).all():
                try:
                    r = BioCReader(source=pubtator.content)
                    r.read()

                    for d_idx, document in enumerate(r.collection.documents):
                        for p_idx, passage in enumerate(document.passages):
                            for annotation in r.collection.documents[d_idx].passages[p_idx].annotations:
                                types_arr.append( annotation.infons['type'] )

                except Exception as e:
                    '''
                    print '%' in pubtator.content
                    for sec in pubtator.document.available_sections():
                        print '%' in sec.text
                    print ' - - - '
                    '''
                    errors = errors + 1

            print 'Errors:', errors
            print Counter(types_arr)
