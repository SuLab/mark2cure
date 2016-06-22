from django.db import models
from django.contrib.auth.models import User

from nltk.tokenize import WhitespaceTokenizer
from mark2cure.common.bioc import BioCReader, BioCWriter, BioCDocument, BioCPassage, BioCAnnotation, BioCLocation

from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes.fields import GenericForeignKey

from ..task.relation.models import Concept, Relation
from ..task.entity_recognition.models import EntityRecognitionAnnotation

import pandas as pd
pd.set_option('display.width', 1000)
import itertools

class Document(models.Model):
    document_id = models.IntegerField(blank=True)
    title = models.TextField(blank=False)
    authors = models.TextField(blank=False)

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)
    source = models.CharField(max_length=200, blank=True)

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

    def as_pubtator_annotation_df(self):
        # If the document has 3 solid annotations
        # "GNormPlus"
        # "DNorm"
        # "tmChem"
        df_columns = ('uid', 'source', 'ann_type', 'text', 'offset', 'location')

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
                for passage in bioc_document.passages:

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

                        #print infons.keys()
                        #print infons
                        #print uid_type, uid, '('+str(annotation_type)+')'
                        #print ' - '*40

                        pubtator_arr.append({
                            'uid': uid,
                            'source': uid_type,

                            'ann_type': annotation_type,
                            'text': str(annotation.text),

                            'offset': int(passage.offset),
                            'location': str(annotation.locations[0])
                        })

                pubtator_dfs.append( pd.DataFrame(pubtator_arr, columns=df_columns) )

        if len(pubtator_dfs):
            return pd.concat(pubtator_dfs)
        else:
            return pd.DataFrame([], columns=df_columns)


    def as_bioc_with_user_annotations(self, request=None):
        '''
            Get the user annotations for a Document
        '''
        document = self.as_bioc()
        approved_types = ['disease', 'gene_protein', 'drug']

        passage_offset = 0
        content_type_id = str(ContentType.objects.get_for_model(EntityRecognitionAnnotation.objects.all().first()).id)
        for section in self.available_sections():
            passage = section.as_bioc(passage_offset)

            for ann in EntityRecognitionAnnotation.objects.annotations_for_section_pk(section.pk, content_type_id):
                annotation = BioCAnnotation()
                annotation.id = str(ann.pk)
                annotation.put_infon('user', str(ann.user_id))

                # (TODO) Map type strings back to 0,1,2
                annotation.put_infon('type', str(approved_types.index(ann.type)))
                annotation.put_infon('type_name', str(ann.type))

                location = BioCLocation()
                location.offset = str(passage_offset + ann.start)
                location.length = str(len(ann.text))
                annotation.add_location(location)

                annotation.text = ann.text
                passage.add_annotation(annotation)

            # (WARNING) Addresses Github Issue #133 & #183
            # Unclear how pubtator will behave with 3+ section documents
            passage_offset += len(passage.text) + 1
            document.add_passage(passage)

        return document

    def as_bioc_with_pubtator_annotations(self, request=None):
        '''
            This is a function that merges the 3 different pubtator
            reponses into 1 main file. It performances selective
            ordering and precedence for some annotations types / instances
        '''
        approved_types = ['Disease', 'Gene', 'Chemical']
        self.init_pubtator()
        reader = self.as_writer(request)

        pub_query_set = Pubtator.objects.filter(
            document=self,
            session_id='',
            content__isnull=False)

        # Load up our various pubtator responses
        pub_readers = []
        for pubtator in pub_query_set.all():
            r = BioCReader(source=pubtator.content)
            r.read()
            pub_readers.append(r)

        for d_idx, document in enumerate(reader.collection.documents):
            for p_idx, passage in enumerate(document.passages):
                # For each passage in each document in the collection
                # add the appropriate annotation
                for p in pub_readers:

                    for annotation in p.collection.documents[d_idx].passages[p_idx].annotations:
                        ann_type = annotation.infons['type']
                        infons = annotation.infons

                        if ann_type in approved_types:

                            uid_type = None
                            uid = None
                            for key in infons.keys():
                                if key != 'type':
                                    uid_type = key
                                    uid = infons.get(uid_type, None)

                            annotation.clear_infons()
                            annotation.put_infon('type', str(approved_types.index(ann_type)))
                            annotation.put_infon('user', 'pubtator')
                            annotation.put_infon('uid', str(uid))
                            reader.collection.documents[d_idx].passages[p_idx].add_annotation(annotation)

                # Remove the shorter annotation if they're multiple
                # at the same start position
                anns = reader.collection.documents[d_idx].passages[p_idx].annotations
                ann_offsets = [x.locations[0].offset for x in anns]

                import collections
                # For each of the offset positions where there are multiple annotations
                for offset in [x for x, y in collections.Counter(ann_offsets).items() if y > 1]:

                    conflicting_anns = [x for x in anns if x.locations[0].offset == offset]
                    longest_ann = max(conflicting_anns, key=lambda a: int(a.locations[0].length))

                    for ann in conflicting_anns:
                        if ann is not longest_ann:
                            reader.collection.documents[d_idx].passages[p_idx].remove_annotation(ann)

                # Remove any annoations that overlap, prefer selection for longest
                anns = reader.collection.documents[d_idx].passages[p_idx].annotations
                for needle_ann in anns:
                    needle_ann_offset = int(needle_ann.locations[0].offset)
                    needle_ann_length = int(needle_ann.locations[0].length)

                    for stack_ann in anns:
                        stack_ann_offset = int(stack_ann.locations[0].offset)
                        stack_ann_length = int(stack_ann.locations[0].length)

                        if needle_ann_offset >= stack_ann_offset and needle_ann_length < stack_ann_length:
                            try:
                                reader.collection.documents[d_idx].passages[p_idx].remove_annotation(needle_ann)
                            except:
                                pass

        return reader.collection.documents[0]


    def as_writer(self, request=None):
        from mark2cure.common.formatter import bioc_writer
        writer = bioc_writer(request)
        document = self.as_bioc_with_passages()
        writer.collection.add_document(document)
        return writer

    def as_bioc_with_passages(self):
        document = self.as_bioc()

        passage_offset = 0
        for section in self.available_sections():
            passage = section.as_bioc(passage_offset)
            passage_offset += len(passage.text)
            document.add_passage(passage)
        return document

    def as_bioc(self):
        document = BioCDocument()
        document.id = str(self.document_id)
        document.put_infon('id', str(self.pk))
        document.put_infon('source', str(self.source))
        return document

    # Helpers for Talk Page
    def annotations(self):
        return EntityRecognitionAnnotation.objects.annotations_for_document_pk(self.pk)

    def contributors(self):
        user_ids = list(set(View.objects.filter(section__document=self, completed=True, task_type='cr').values_list('user', flat=True)))
        return user_ids



    class Meta:
        ordering = ('-created',)
        get_latest_by = 'updated'


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
        # (TODO) This is not supported by MySQL but would help prevent dups in this table
        # unique_together = ['kind', 'type', 'text', 'start', 'view']
