from django.core.management.base import BaseCommand

from ..document.models import Pubtator
from ..common.bioc import BioCReader

from optparse import make_option
from collections import Counter


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
                                types_arr.append(annotation.infons['type'])

                except Exception:
                    '''
                    print '%' in pubtator.content
                    for sec in pubtator.document.available_sections():
                        print '%' in sec.text
                    print ' - - - '
                    '''
                    errors = errors + 1

            print 'Errors:', errors
            print Counter(types_arr)
