from django.db import models
from django.contrib.auth.models import User

from nltk.tokenize import WhitespaceTokenizer
from mark2cure.common.bioc import BioCReader, BioCWriter, BioCPassage
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes.fields import GenericForeignKey

from .managers import DocumentManager
from ..task.entity_recognition.models import EntityRecognitionAnnotation
# from ..task.relation.models import RelationAnnotation

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

    def init_pubtator(self):
        if self.available_sections().exists() and Pubtator.objects.filter(document=self).count() < 3:
            for api_ann in ['tmChem', 'DNorm', 'GNormPlus']:
                Pubtator.objects.get_or_create(document=self, kind=api_ann)

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
        # Quick check using cache
        if Pubtator.objects.filter(document=self, validate_cache=True).count() == 3:
            return True

        # All responses that are not waiting to fetch conent b/c they already have it
        pubtators = Pubtator.objects.filter(
            document=self,
            session_id='',
            content__isnull=False).all()

        # Check if each type validates
        if pubtators.count() != 3:
            return False

        for pubtator in pubtators:
            if not pubtator.valid():
                return False

        return True

    DF_COLUMNS = ('uid', 'source', 'user_id',
                  'ann_type', 'text',
                  'section_id', 'section_offset', 'offset_relative',
                  'start_position', 'length')
    APPROVED_TYPES = ['disease', 'gene_protein', 'drug']

    def create_er_df_row(self,
                      uid, source='db', user_id=None,
                      ann_type='', text='',
                      section_id=0, section_offset=0, offset_relative=True,
                      start_position=0, length=0):
        '''
            When offset_relative is False:
                start position is relative to the entire document and not the
                section it's contained within

            user_id can be None

            (TODO) Ann Types needs to be normalized
        '''
        ann_type = ann_type.lower()
        return {
            'uid': str(uid), 'source': str(source), 'user_id': int(user_id) if user_id else None,
            'ann_type': str(ann_type), 'text': str(text),
            'section_id': int(section_id), 'section_offset': int(section_offset), 'offset_relative': bool(offset_relative),
            'start_position': int(start_position), 'length': int(length)
        }

    def as_er_df_with_pubtator_annotations(self):
        '''
            This is a function that merges the 3 different pubtator
            reponses into 1 main file. It performances selective
            ordering and precedence for some annotations types / instances
        '''
        # If the document has 3 solid annotations
        # "GNormPlus"
        # "DNorm"
        # "tmChem"

        section_ids = self.section_set.values_list('pk', flat=True)

        pubtator_dfs = []
        if self.valid_pubtator():

            pubtators = Pubtator.objects.filter(
                document=self,
                session_id='',
                content__isnull=False).all()

            for pubtator in pubtators:
                r = BioCReader(source=pubtator.content)
                r.read()

                pubtator_arr = []
                bioc_document = r.collection.documents[0]
                for p_idx, passage in enumerate(bioc_document.passages):

                    for annotation in passage.annotations:
                        infons = annotation.infons

                        annotation_type = None
                        uid_type = None
                        uid = None

                        for key in infons.keys():
                            if key == 'type':
                                annotation_type = infons.get(key, None)
                            else:
                                uid_type = key
                                uid = infons.get(uid_type, None)

                        start, length = str(annotation.locations[0]).split(':')
                        pubtator_arr.append(self.create_er_df_row(
                            uid=uid, source=uid_type, user_id=None,
                            text=annotation.text, ann_type=annotation_type,
                            section_id=section_ids[p_idx], section_offset=passage.offset, offset_relative=False,
                            start_position=start, length=length))

                pubtator_dfs.append(pd.DataFrame(pubtator_arr, columns=Document.DF_COLUMNS))

        if len(pubtator_dfs):
            return pd.concat(pubtator_dfs)
        else:
            return pd.DataFrame([], columns=Document.DF_COLUMNS)

    def as_er_df_with_user_annotations(self, user=None):
        '''
            Returns back a Pandas Dataframe with the Entity Recognition annotations
            submitted by all users for this document.

            If a user is passed in, the returning dataframe only contains annotations
            by that individual user.
        '''

        df_arr = []

        content_type_id = str(ContentType.objects.get_for_model(
            EntityRecognitionAnnotation.objects.first()).id)

        if user:
            # (id, type, text, start, created, section_id, user_id)
            er_ann_query_set = EntityRecognitionAnnotation.objects.\
                annotations_for_document_pk_and_user(self.pk, user.pk, content_type_id)
        else:
            # (id, type, text, start, created, section_id, user_id)
            er_ann_query_set = EntityRecognitionAnnotation.objects.\
                annotations_for_document_pk(self.pk, content_type_id)

        # Use the BioC pubtator file for the offset values
        offset_dict = {}
        writer = self.as_writer()
        for passage in writer.collection.documents[0].passages:
            offset_dict[int(passage.infons.get('id'))] = passage.offset

        for er_ann in er_ann_query_set:

            df_arr.append(self.create_er_df_row(
                uid=er_ann.id, source='db', user_id=er_ann.user_id,
                text=er_ann.text, ann_type=er_ann.type,
                section_id=er_ann.section_id, section_offset=offset_dict[er_ann.section_id], offset_relative=True,
                start_position=er_ann.start, length=len(er_ann.text)))

        return pd.DataFrame(df_arr, columns=Document.DF_COLUMNS)

    def as_writer(self):
        '''
            Return a blank BioC Writer that is based off the pubtator content.

            Problems: This requires every document to have at least 1 pubtator model
            Pros: This prevents us from generating our own BioC file which may
            have inconsistencies
        '''
        pubtator = Pubtator.objects.filter(
            document=self,
            session_id='',
            content__isnull=False).first()

        if not pubtator:
            return False

        r = BioCReader(source=pubtator.content)
        r.read()

        section_ids = self.section_set.values_list('pk', flat=True)

        for doc in r.collection.documents:
            for idx, passage in enumerate(doc.passages):
                passage.clear_annotations()

                passage.put_infon('section', ['title', 'paragraph'][idx])
                passage.put_infon('id', str(section_ids[idx]))

        bioc_writer = BioCWriter()
        bioc_writer.collection = r.collection

        return bioc_writer

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
    document = models.ForeignKey(Document)

    kind = models.CharField(max_length=200, blank=True)
    session_id = models.CharField(max_length=200, blank=True)
    content = models.TextField(blank=True, null=True)

    request_count = models.IntegerField(default=0)
    validate_cache = models.BooleanField(default=False)

    # Updated is also trigged during doc.valid_pubtator
    # so it's associated with when last polled or
    # since last time it's been known to validate
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    class Meta:
        app_label = 'document'

    def __unicode__(self):
        return 'pubtator'

    def valid(self):
        # (TODO) This may return 2 different "types" check on
        # implications of this discrepancy
        if self.validate_cache:
            return True

        if self.session_id != '':
            return False

        if self.content is None:
            return False

        try:
            r = BioCReader(source=self.content)
            r.read()
            return r
        except Exception:
            # If one of them doesn't validate leave
            return False

    def as_writer(self, request=None):
        r = BioCReader(source=self.content)
        r.read()

        bioc_writer = BioCWriter()
        bioc_writer.collection = r.collection

        return bioc_writer

    def count_annotations(self):
        if self.valid():
            count = 0
            reader = BioCReader(source=self.content)
            reader.read()
            for doc in reader.collection.documents:
                for passage in doc.passages:
                    count += len(passage.annotations)
            return count

        else:
            return 0


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

    def as_bioc(self, offset, content=True):
        passage = BioCPassage()
        passage.put_infon('type', 'paragraph')
        passage.put_infon('section', self.get_kind_display().lower())
        passage.put_infon('id', str(self.pk))
        if content:
            passage.text = self.text
        else:
            passage.text = ''
        passage.offset = str(offset)
        return passage

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
        return u'Document #{doc_id}, Section #{sec_id} by {username}'.format(
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
        if self.kind == 'e':
            return '{0} ({1}) [{2}]'.format(self.text, self.start, self.pk)
        else:
            return self.type

    class Meta:
        get_latest_by = 'updated'
        app_label = 'document'

        # (TODO) This is not supported by MySQL but would help prevent dups in this table
        # unique_together = ['kind', 'type', 'text', 'start', 'view']
