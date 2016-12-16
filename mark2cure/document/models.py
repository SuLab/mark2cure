from django.db import models
from django.contrib.auth.models import User
from django.utils import timezone

from nltk.tokenize import WhitespaceTokenizer
from mark2cure.common.bioc import BioCReader
from mark2cure.common.formatter import validate_pubtator
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes.fields import GenericForeignKey

from .managers import DocumentManager
from ..task.entity_recognition.models import EntityRecognitionAnnotation
from librabbitmq import ConnectionError

import pandas as pd
pd.set_option('display.width', 1000)


class Document(models.Model):
    document_id = models.IntegerField(blank=True)
    title = models.TextField(blank=False)
    authors = models.TextField(blank=False)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)
    source = models.CharField(max_length=200, blank=True)

    objects = DocumentManager()

    def __unicode__(self):
        return self.title

    def available_sections(self):
        return self.section_set.exclude(kind='o').all()

    def count_available_sections(self):
        return self.section_set.exclude(kind='o').count()

    def run_pubtator(self):
        """ Ensure a Document has Pubtator entries
            1) Create Pubtator entries
            2) Start checking for responses if pending
        """
        kind_arr = ['tmChem', 'DNorm', 'GNormPlus']

        # If the Document has excess Pubtators, cleanup
        if self.pubtators.count() > 3:
            for kind in kind_arr:
                if self.pubtators.filter(kind=kind).count() > 1:
                    # Picks them off 1 run at a time (oldest first) until no duplicates remain
                    self.pubtators.filter(kind=kind).first().delete()

        # If the Document is missing any Pubtator "kinds", it should get them
        for missing_kind in list(set(kind_arr) - set(self.pubtators.values_list('kind', flat=True))):
            Pubtator.objects.get_or_create(document=self, kind=missing_kind)

        for pubtator in self.pubtators.filter(content__isnull=True).all():
            try:
                pending_request = pubtator.requests.get(status=PubtatorRequest.UNFULLFILLED)

                # If the current request is never going to finish, flag it and start over
                if pending_request and (timezone.now() - pending_request.updated).days >= 1:
                    pending_request.status = PubtatorRequest.EXPIRED
                    pending_request.save()
                    pubtator.submit()

            except PubtatorRequest.DoesNotExist:
                # Start the PubtatorRequest for the first time
                pubtator.submit()

            except PubtatorRequest.MultipleObjectsReturned:
                # If multiple arise, delete all and restart the PubtatorRequest attempts
                pubtator.requests.filter(status=PubtatorRequest.UNFULLFILLED).all().delete()
                pubtator.submit()

    def update_padding(self):
        from mark2cure.common.formatter import pad_split
        changed = False

        for section in self.available_sections():
            padded = ' '.join(pad_split(section.text))
            if section.text != padded:
                # If a change was identified:
                # 1) Resubmit it to pubtator
                # 2) Remove any submissions for this doc OR flag their annotations
                section.text = padded
                section.save()
                changed = True

        return changed

    def valid_pubtator(self):
        """ Returns a boolean if the Document instance has a complet Pubtator
            coverage for all kinds (tmChem, DNorm, GNormPlus)
        """
        # All responses that are not waiting to fetch conent b/c they already have it
        pubtators = self.pubtators.all()

        # Check if each type validates
        if pubtators.count() != 3:
            return False

        for pubtator in pubtators:
            if not validate_pubtator(pubtator.content, pubtator.document):
                return False

        return True

    # Helpers for Talk Page
    def annotations(self):
        return EntityRecognitionAnnotation.objects.annotations_for_document_pk(self.pk)

    def contributors(self):
        user_ids = list(set(View.objects.filter(section__document=self, completed=True, task_type='cr').values_list('user', flat=True)))
        return user_ids

    class Meta:
        ordering = ('-created',)
        get_latest_by = 'updated'
        app_label = 'document'


class Pubtator(models.Model):
    document = models.ForeignKey(Document, related_name='pubtators')
    kind = models.CharField(max_length=200)
    content = models.TextField(blank=True, null=True)

    class Meta:
        app_label = 'document'

    def __unicode__(self):
        return '{0} for PMID: {1} ({2})'.format(self.kind, self.document.document_id, 'Valid' if self.is_valid() else 'Invalid')

    def is_valid(self):
        """ Alias for validate_pubtator
            Returns a boolean for if the pubtator is a valid state
        """
        return validate_pubtator(self.content, self.document)

    def count_annotations(self):
        """ Returns an int count of all types of ER annotations in the Pubtator instance
            If none are found or the document is invalid, return 0
        """
        try:
            r = BioCReader(source=self.content)
            r.read()
            return sum([len(passage.annotations) for passage in r.collection.documents[0].passages])
        except Exception:
            return 0

    def submit(self):
        from .tasks import submit_pubtator
        try:
            submit_pubtator.apply_async(
                args=[self.pk],
                queue='mark2cure_tasks')
        except ConnectionError:
            submit_pubtator(self.pk)


class PubtatorRequest(models.Model):
    """ Pending jobs that have been submitted to Pubtator and are
        awaiting completion
    """
    pubtator = models.ForeignKey(Pubtator, related_name='requests')

    UNFULLFILLED = 0
    FULLFILLED = 1
    FAILED = 2
    EXPIRED = 3
    STATUS_CHOICES = (
        (UNFULLFILLED, 'Unfullfilled'),
        (FULLFILLED, 'Fullfilled'),
        (FAILED, 'Failed'),
        (EXPIRED, 'Expired')
    )
    status = models.IntegerField(default=UNFULLFILLED, choices=STATUS_CHOICES)
    session_id = models.CharField(max_length=19)
    # The Number of times we've checked on the session_id
    request_count = models.IntegerField(default=0)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    def __unicode__(self):
        return '{0} for Pubtator ({1})'.format(self.session_id, self.pubtator)

    def check_status(self):
        from .tasks import check_pubtator
        try:
            check_pubtator.apply_async(
                args=[self.pk],
                queue='mark2cure_tasks')
        except ConnectionError:
            check_pubtator(self.pk)


class Section(models.Model):
    SECTION_KIND_CHOICE = (
        ('o', 'Overview'),
        ('t', 'Title'),
        ('a', 'Abstract'),
        ('p', 'Paragraph'),
        ('f', 'Figure'),
    )
    kind = models.CharField(max_length=1, choices=SECTION_KIND_CHOICE)
    text = models.TextField(blank=True)
    source = models.ImageField(blank=True, upload_to='media/images/', default='images/figure.jpg')

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    document = models.ForeignKey(Document)

    class Meta:
        app_label = 'document'

    def resultwords(self, user_view, gm_view):
        # Gather words and positions from the text
        words_index = WhitespaceTokenizer().span_tokenize(self.text)
        words_text = WhitespaceTokenizer().tokenize(self.text)
        words = zip(words_index, words_text)

        # Add counters for concensus count and personal annotation
        # ((start, stop), 'Word string itself', Intensity, GM Ann ID, User Ann ID, Did user annotate)
        words = [w + (0, None, None, False,) for w in words]

        # Gather other annotations from GM and users for this section
        gm_anns = Annotation.objects.filter(view=gm_view).values_list('pk', 'start', 'text')

        # Build the running counter of times a word was annotated
        for gm_pk, start, text in gm_anns:
            length = len(text)

            for idx, word in enumerate(words):
                word_start = word[0][0]
                counter = word[2]
                if word_start >= start and word_start <= start + length:
                    counter += 1
                    words[idx] = (word[0], word[1], counter, gm_pk, word[3], word[4])

        user_anns = Annotation.objects.filter(view=user_view).values_list('pk', 'start', 'text')

        # Build the running counter of times a word was annotated
        for user_pk, start, text in user_anns:
            length = len(text)
            for idx, word in enumerate(words):
                word_start = word[0][0]
                if word_start >= start and word_start <= start + length:
                    words[idx] = (word[0], word[1], word[2], word[3], user_pk, True)

        return words

    def update_view(self, user, task_type, completed=False):
        view = View.objects.filter(user=user, task_type=task_type, section=self).latest()
        view.completed = completed
        view.save()
        return view

    def __unicode__(self):
        if self.kind == 'o':
            return u'[Overview] {0}'.format(self.document.title)
        else:
            return self.text

TASK_TYPE_CHOICE = (
    ('cr', 'Concept Recognition'),
    ('cn', 'Concept Normalization'),
    ('rv', 'Relationship Verification'),
    ('ri', 'Relationship Identification'),
    ('rc', 'Relationship Correction'),
)


class View(models.Model):
    '''
        When completing tasks not on a Section level. Work is associated
        with the FIRST section available, regardless of the Section type
    '''
    task_type = models.CharField(max_length=3, choices=TASK_TYPE_CHOICE, blank=True, default='cr')
    completed = models.BooleanField(default=False, blank=True)
    opponent = models.ForeignKey('self', blank=True, null=True)

    section = models.ForeignKey(Section)
    user = models.ForeignKey(User)

    class Meta:
        get_latest_by = 'pk'
        app_label = 'document'

    def __unicode__(self):
        return u'{pk}, Document #{doc_id}, Section #{sec_id} by {username}'.format(
            pk=self.pk,
            doc_id=self.section.document.pk,
            sec_id=self.section.pk,
            username=self.user.username)


class Annotation(models.Model):
    ANNOTATION_KIND_CHOICE = (
        ('e', 'Entity Recognition'),
        ('r', 'Relation'),
    )
    kind = models.CharField(max_length=1, choices=ANNOTATION_KIND_CHOICE, blank=False, default='e')

    created = models.DateTimeField(auto_now_add=True)

    content_type = models.ForeignKey(ContentType, blank=True, null=True)
    object_id = models.IntegerField(blank=True, null=True)
    metadata = GenericForeignKey('content_type', 'object_id')

    view = models.ForeignKey(View, blank=True, null=True)

    def __unicode__(self):
        if self.kind == 'r':
            return 'Relationship Ann'
        elif self.kind == 'e':
            return '{0} ({1}) [{2}]'.format(self.text, self.start, self.pk)
        else:
            return 'Annotation {0}'.format(self.pk)

    class Meta:
        get_latest_by = 'updated'
        app_label = 'document'
