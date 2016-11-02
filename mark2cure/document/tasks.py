from __future__ import absolute_import

from django.conf import settings

from mark2cure.common.models import Group
from mark2cure.document.models import Document, Pubtator, PubtatorRequest, Section
from mark2cure.common.formatter import pad_split

from Bio import Entrez, Medline
from ..common import celery_app as app

from celery.exceptions import SoftTimeLimitExceeded

import requests
import logging
import re
logger = logging.getLogger(__name__)


@app.task(bind=True, ignore_result=True,
          max_retries=0,
          acks_late=True, track_started=True,
          expires=300)
def check_corpus_health(self):
    """
        Task to run every 10 minutes

        1) Check on document health
            - Content is present
            - Padding and Pubtator on new content

        2) Check on the pubtator health
            - Fetch pending sessions
            - Make sure content pubtators are correctly assigned
    """
    for document in Document.objects.all():
        # Update any documents that don't have a Title or Abstract
        if document.available_sections().count() < 2:
            get_pubmed_document.apply_async(
                args=[document.document_id],
                queue='mark2cure_tasks')

            # Update any newly enforced padding rules
            # If the document doesn't pass the validator
            # delete all existing content and retry
            if not document.valid_pubtator() or document.update_padding():
                document.init_pubtator()

    if not self.request.called_directly:
        return True


@app.task(bind=True, ignore_result=True,
          max_retries=0, rate_limit='2/s', soft_time_limit=15,
          acks_late=True, track_started=True,
          expires=60)
def get_pubmed_document(self, pubmed_ids, source='pubmed', include_pubtator=True, group_pk=None):
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

    if not self.request.called_directly:
        return True


@app.task(bind=True, ignore_result=True,
          max_retries=0,
          acks_late=True, track_started=True,
          expires=180)
def maintain_pubtator_requests(self):
    """A routine job that continually checks for pending Pubtator Requests
        and will resubmit Pubtators that haven't been updated in a "long time"
    """

    # Try to fetch all the pending pubtator requests
    for pubtator_request in PubtatorRequest.objects.filter(fulfilled=False).all():
        check_pubtator.apply_async(
            args=[pubtator_request.pk],
            queue='mark2cure_tasks')

    if not self.request.called_directly:
        return True


@app.task(bind=True, ignore_result=True,
          max_retries=0, rate_limit='2/s', soft_time_limit=15,
          acks_late=True, track_started=True,
          expires=None)
def submit_pubtator(self, pubtator_pk):
    """Takes an existing Pubtator instance and submits a processing request
    """
    pubtator = Pubtator.objects.get(pk=pubtator_pk)

    # Make response to post job to pubtator
    payload = {'content-type': 'text/xml'}
    writer = Document.objects.as_writer(documents=[pubtator.document])
    data = str(writer)
    url = 'http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/{api_ann}/Submit/'.format(api_ann=pubtator.kind)

    try:
        response = requests.post(url, data=data, params=payload)
    except Exception as e:
        raise e

    session_id = re.findall(r'\d{4}-\d{4}-\d{4}-\d{4}', response.url)[0]
    PubtatorRequest.objects.get_or_create(
        pubtator=pubtator,
        session_id=session_id)

    if not self.request.called_directly:
        return True


@app.task(bind=True, ignore_result=True,
          max_retries=0, rate_limit='2/s', soft_time_limit=15,
          acks_late=True, track_started=True,
          expires=None)
def check_pubtator(self, pubtator_request_pk):
    """Takes a Pubtator Request and checks for the status from the server
    """
    pubtator_request = PubtatorRequest.objects.get(pk=pubtator_request_pk)
    pubtator = pubtator_request.pubtator

    # Build the writer required to fetch content from previous session
    # Unclear why Pubtator wants this posted twice
    writer = Document.objects.as_writer(documents=[pubtator.document])

    try:
        results = requests.post('http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/{session_id}/Receive/'.format(session_id=pubtator_request.session_id),
                                data=str(writer),
                                params={'content-type': 'text/xml'})
        pubtator_request.request_count = pubtator_request.request_count + 1
        pubtator_request.save()
    except SoftTimeLimitExceeded:
        return False

    except:
        return False

    if results.content != 'Not yet':
        pubtator.content = results.text
        pubtator.save()

        pubtator_request.fulfilled = True
        pubtator_request.save()

    if not self.request.called_directly:
        return True

