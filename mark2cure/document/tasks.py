from django.conf import settings

from mark2cure.common.models import Group
from mark2cure.document.models import Document, Pubtator, Section
from mark2cure.common.formatter import pad_split

from Bio import Entrez, Medline
from celery import task

import requests


@task()
def check_system_uptime():
    logger.info('Celery on rabbitmq-server ran.')


@task()
def check_corpus_health(group_pk=1):
    group = Group.objects.get(pk=group_pk)
    documents = Document.objects.filter(task__group=group).all()

    for document in documents:
        # Update any documents that don't have a Title or Abstract
        if document.available_sections().count() < 2:
            print '  > Downloading Sections: '
            get_pubmed_document(document.document_id)

        # Update any newly enforced padding rules
        # If the document doesn't pass the validator
        # delete all existing content and retry
        valid_pubtator_responses = document.valid_pubtator()
        update_padding = document.update_padding()

        print '  > Valid Pub: ', valid_pubtator_responses
        print '  > Updated Padding: ', update_padding

        if not valid_pubtator_responses or update_padding:
            # Delete only the Pubtator models that are not currently
            # being fetched
            # Pubtator.objects.filter(
            #    document=document,
            #    session_id='',
            #    content__isnull=True).all().delete()

            document.init_pubtator()
            print ' > New Pubtator requests for Doc #', document.pk

        print 'Done checking Doc #', document.pk


def check_pubtator_health():
    # Try to fetch all the pending pubtator requests
    for pubtator in Pubtator.objects.exclude(session_id='').all():
        get_pubtator_response(pubtator.pk)

    # For all Pubtator models with content
    # ensure it validates and cleanup the session_id and content
    # if it does or doesn't
    Pubtator.objects.correct_parent_relation()


@task()
def get_pubtator_response(pk):
    pubtator = Pubtator.objects.get(pk=pk)

    if pubtator.session_id:
        # Body required to fetch content from previous session
        payload = {'content-type': 'text/xml'}
        writer = pubtator.document.as_writer()
        data = str(writer)
        url = 'http://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/{session_id}/Receive/'.format(
            session_id=pubtator.session_id)

        results = requests.post(url, data=data, params=payload)
        pubtator.request_count = pubtator.request_count + 1

        if results.content != 'Not yet':
            pubtator.session_id = ''
            pubtator.content = results.text

        pubtator.save()


@task()
def get_pubmed_document(pubmed_ids, source='pubmed', include_pubtator=True):
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
            # if Document.objects.pubmed_count(record.get('PMID')) is 0:
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


@task
def get_pubmed_documents(terms=settings.ENTREZ_TERMS):
    Entrez.email = settings.ENTREZ_EMAIL

    for term in terms:
        h = Entrez.esearch(db='pubmed', retmax=settings.ENTREZ_MAX_COUNT, term=term)
        result = Entrez.read(h)
        ids = result['IdList']
        get_pubmed_document(ids)

