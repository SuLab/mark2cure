from django.db import models
from django.contrib.auth.models import User

from mark2cure.document.managers import DocumentManager, PubtatorManager

from nltk.tokenize import WhitespaceTokenizer
from mark2cure.common.bioc import BioCReader, BioCDocument, BioCPassage, BioCAnnotation, BioCLocation
from django.forms.models import model_to_dict


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
        # All responses that are not waiting to fetch conent b/c they already have it
        pub_query_set = Pubtator.objects.filter(
            document=self,
            session_id='',
            content__isnull=False)

        for section in self.available_sections():
            for needle in ['<', '>']:
                if needle in section.text:
                    return False

        # The Docment doesn't have a response for each type
        # (TODO) also cases grater than 3
        if pub_query_set.count() != 3:
            return False

        # Check if each type validates, if so save
        for pubtator in pub_query_set.all():

            p_valid = pubtator.valid()
            if p_valid:
                pubtator.document = Document.objects.get(document_id=p_valid.collection.documents[0].id)
                # This updates the created field, so we know the last time it was checked
                pubtator.save()
            else:
                return False
        return True

    def as_bioc_with_user_annotations(self, request=None):
        '''
            Get the user annotations for a Document
        '''
        document = self.as_bioc()
        approved_types = ['disease', 'gene_protein', 'drug']

        passage_offset = 0
        for section in self.available_sections():
            passage = section.as_bioc(passage_offset)

            for ann in Annotation.objects.filter(view__section=section).values('pk', 'start', 'text', 'type', 'view__user__username', 'view__user__pk').all():
                annotation = BioCAnnotation()
                annotation.id = str(ann.get('pk'))
                annotation.put_infon('user', str(ann.get('view__user__pk')))
                annotation.put_infon('user_name', str(ann.get('view__user__username')))

                # (TODO) Map type strings back to 0,1,2
                annotation.put_infon('type', str(approved_types.index(ann.get('type'))))
                annotation.put_infon('type_name', str(ann.get('type')))

                location = BioCLocation()
                location.offset = str(passage_offset + ann.get('start'))
                location.length = str(len(ann.get('text')))
                annotation.add_location(location)

                annotation.text = ann.get('text')
                passage.add_annotation(annotation)

            passage_offset += len(passage.text)
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

                        if ann_type in approved_types:
                            annotation.clear_infons()
                            annotation.put_infon('type', str(approved_types.index(ann_type)))
                            annotation.put_infon('user', 'pubtator')
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
        return Annotation.objects.filter(view__section__document=self)

    def contributors(self):
        user_ids = list(set(View.objects.filter(section__document=self, completed=True).values_list('user', flat=True)))
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

    objects = PubtatorManager()

    def __unicode__(self):
        return 'pubtator'

    def valid(self):
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

    def count_annotations(self):
        reader = self.valid()
        count = 0

        if reader:
            for doc in reader.collection.documents:
                for passage in doc.passages:
                    count += len(passage.annotations)

        return count


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

    def as_bioc(self, offset):
        passage = BioCPassage()
        passage.put_infon('type', 'paragraph')
        passage.put_infon('section', self.get_kind_display().lower())
        passage.put_infon('id', str(self.pk))
        passage.text = self.text
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
        ('e', 'Entities'),
        ('a', 'Attributes'),
        ('r', 'Relations'),
        ('t', 'Triggers'),
    )
    kind = models.CharField(max_length=1, choices=ANNOTATION_KIND_CHOICE, blank=False, default='e')

    # Disease, Gene, Protein, et cetera...
    ANNOTATION_TYPE_CHOICE = (
        'disease',
        'gene_protein',
        'drug',
    )
    type = models.CharField(max_length=40, blank=True, null=True, default='disease')

    text = models.TextField(blank=True, null=True)

    # (WARNING) Different than BioC
    # This is always the start position relative
    # to the section, not the entire document
    start = models.IntegerField(blank=True, null=True)

    created = models.DateTimeField(auto_now_add=True)

    view = models.ForeignKey(View)

    def is_exact_match(self, comparing_annotation):
        required_matching_keys = ['kind', 'start', 'text', 'type']
        self_d = model_to_dict(self)
        compare_d = model_to_dict(comparing_annotation)
        return all([True if self_d[k] == compare_d[k] else False for k in required_matching_keys])

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

