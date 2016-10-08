from __future__ import absolute_import

from django.conf import settings

from mark2cure.common.models import Group
from mark2cure.document.models import Document, Pubtator, Section
from mark2cure.common.formatter import pad_split

from Bio import Entrez, Medline
from ..common import celery_app as app

import requests
import logging
logger = logging.getLogger(__name__)


@app.task(bind=True, ignore_result=False,
          max_retries=1,
          acks_late=True, track_started=True,
          expires=None)
def check_corpus_health(self):
    '''
        Task to run every 10 minutes

        1) Check on document health
            - Content is present
            - Padding and Pubtator on new content

        2) Check on the pubtator health
            - Fetch pending sessions
            - Make sure content pubtators are correctly assigned

    '''
    for document in Document.objects.all():
        # Update any documents that don't have a Title or Abstract
        if document.available_sections().count() < 2:
            get_pubmed_document(document.document_id)

            # Update any newly enforced padding rules
            # If the document doesn't pass the validator
            # delete all existing content and retry
            valid_pubtator_responses = document.valid_pubtator()
            update_padding = document.update_padding()

            if not valid_pubtator_responses or update_padding:
                document.init_pubtator()

    check_pubtator_health()


def check_pubtator_health():
    # Set Validate Cache to False for all to perform
    # an entire, clean sweep of new checks
    Pubtator.objects.all().update(validate_cache=False)

    # Try to fetch all the pending pubtator requests
    for pubtator in Pubtator.objects.exclude(session_id='').all():
        get_pubtator_response(pubtator.pk)

    # For all Pubtator models with content ensure it validates and cleanup the session_id and content
    for pubtator in Pubtator.objects.filter(content__isnull=False).all():
        # (TODO) Do robust checks for Pubtator object valid status
        p_valid = pubtator.valid()

        if p_valid:
            # Association with the correct document
            pubtator.document = Document.objects.get(document_id=p_valid.collection.documents[0].id)

            # Prevents subsequent API calls
            pubtator.session_id = ''

            # Make Pubtator.validate() faster
            pubtator.validate_cache = True

        else:
            pubtator.content = None

        # Do this just so the first time valid_pubtator
        # actually runs we know it's fresh'
        pubtator.save()


@app.task(bind=True, ignore_result=False,
          max_retries=1, rate_limit='2/s', soft_time_limit=15,
          acks_late=True, track_started=True,
          expires=None)
def get_pubtator_response(pk):
    pubtator = Pubtator.objects.get(pk=pk)

    if pubtator.session_id:
        # Body required to fetch content from previous session
        payload = {'content-type': 'text/xml'}
        writer = Document.objects.as_writer(documents=[pubtator.document])
        data = str(writer)
        url = 'http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/{session_id}/Receive/'.format(
            session_id=pubtator.session_id)

        results = requests.post(url, data=data, params=payload)
        pubtator.request_count = pubtator.request_count + 1

        if results.content != 'Not yet':
            pubtator.session_id = ''
            pubtator.content = results.text

        pubtator.save()


@app.task(bind=True, ignore_result=False,
          max_retries=1, rate_limit='2/s', soft_time_limit=15,
          acks_late=True, track_started=True,
          expires=None)
def get_pubmed_document(pubmed_ids, source='pubmed', include_pubtator=True, group_pk=None):
    Entrez.email = settings.ENTREZ_EMAIL

    if type(pubmed_ids) == list:
        ids = [str(doc_id) for doc_id in pubmed_ids]
    else:
        ids = [str(pubmed_ids)]

    h = Entrez.efetch(db='pubmed', id=ids, rettype='medline', retmode='text')
    records = Medline.parse(h)

    # Reference to abbreviations: http://www.nlm.nih.gov/bsd/mms/medlineelements.html
    for record in records:
        if record.get('TI') and record.get('AB') and record.get('PMID') and record.get('CRDT'):
            title = ' '.join(pad_split(record.get('TI')))
            abstract = ' '.join(pad_split(record.get('AB')))

            doc, doc_c = Document.objects.get_or_create(document_id=record.get('PMID'))
            doc.title = title
            doc.source = source
            doc.save()

            sec, sec_c = Section.objects.get_or_create(kind='t', document=doc)
            sec.text = title
            sec.save()

            sec, sec_c = Section.objects.get_or_create(kind='a', document=doc)
            sec.text = abstract
            sec.save()

            if include_pubtator:
                doc.init_pubtator()

    if group_pk:
        docs = Document.objects.filter(source=source).all()
        group = Group.objects.get(pk=group_pk)
        group.assign(docs)
