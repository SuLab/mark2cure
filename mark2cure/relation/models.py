'''
Sets up the models for the relation application, modeled after Max's
Document/models.py file for standardization and using Toby's data_model.py

Most of this code is Toby's, but reconfigured to work in the Django app, where
tasks.py performs much of the actions involved in populating the new documents
(pmid) into the database.

'''
from django.db import models
from ..document.models import Document

class Paper(models.Model):
    """
    snip of the Max's document section in the database administration:
    Document id
    Title preview    Sections   Pubtator   Annotations   Completed views    Pending   views  Source
    22206664	The Nogo receptor 2 is a novel substrate of Fbs1 .	2	True	82	12	2	group5

    """
    #def __init__(self, pmid, title, abstract, annotations, gold_relations = []):
    # 1
    pmid = models.TextField(blank=False)
    # 2
    title = models.TextField(blank=False)
    # 3
    abstract = models.TextField(blank=False)
    # 4
    # annotations = sorted(annotations)
    """ Without having to put the "unique chemicals and unique
    diseases" in this model, put into a function to filter, because
    annotations are already in their own, distinct model  """
    def get_unique_diseases(self):
        # return a list of the unique diseases to this Paper
        query = Concept.objects.filter(paper=self).filter(stype='disease').values_list('uid').distinct()
        return query

    def get_unique_chemicals(self):
        # return a list of the unique chemicals to this Paper
        query = Concept.objects.filter(paper=self).filter(stype='chemical').values_list('uid').distinct()
        return query

    # gold_relations = models.TextField(blank=False) # TODO check what false is here
    # assert self._has_correct_annotations() TODO add back
    # 5
    # TODO add chemicals and diseases and other items to model
    # gold_relations = models.TextField(blank=False) # may be empty when not parsing gold # LIST TODO
    # 6 & 7

    # 8 split sentences and generate sentence-bound relations
    # sentences = _split_sentences()
    # 9
    #possible_relations = _classify_relations()

    # Max's model
    """
    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)
    """
    """
    def available_sections(self):
        return self.section_set.all()

    def count_available_sections(self):
        return self.section_set.count()
    """
    # TODO add more methods that include objects of the below classes for relation module

    """
    # TODO fix the contributors (we want to know users that contributed to this and the below CLASSES)
    def contributors(self):
        user_ids = list(set(View.objects.filter(section__document=self, completed=True).values_list('user', flat=True)))
        return user_ids
    """
    # follows max's model
    def __unicode__(self):
        return self.pmid
    """
    # follows Toby's model:
    def __repr__(self):
        return ("<{0}>: PMID {1}. {2} annotations, {3} gold relations\n"
            "{4} unique chemical ids, {5} unique disease ids\n"
            "{6} sentences".format(self.__class__.__name__,
            self.pmid, len(self.annotations), len(self.gold_relations),
            len(self.chemicals), len(self.diseases), len(self.sentences))
        )
    """
    def _has_correct_annotations(self):
        """
        Checks that the paper's annotations match the
        stated positions in the text.
        """
        text = "{0} {1}".format(self.title, self.abstract)
        for annotation in self.annotations:
            assert text[annotation.start:annotation.stop] == annotation.text, (
                "Annotation {0} in PMID {1} does not match the text.".format(annotation, self.pmid))
        return True

    def _get_unique_concepts(self):
        """
        Determines the unique identifiers of chemicals and diseases
        belonging to this paper.

        Ignores any annotations with an identifier of -1.
        """
        res = defaultdict(set)
        for annotation in self.annotations:
            if annotation.uid != "-1":
                res[annotation.stype].add(annotation.uid)
        return (res["chemical"], res["disease"])

    def _get_all_possible_relations(self):
        """
        Returns all possible unique drug-disease combinations
        as a set for this paper.
        """
        return {(chemical_id, disease_id) for chemical_id in self.chemicals for disease_id in self.diseases}

    def _split_sentences(self):
        """
        Splits the abstract up into individual sentences,
        and determines which concept annotations reside
        within each sentence.

        Time complexity:
            O(N + M) where N is the number of sentences,
            and M is the number of annotations.
        """
        all_sentences = [self.title] + split_abstract(self.abstract)

        full_text = "{0} {1}".format(self.title, self.abstract)

        sent_idx = 0  # starting index of current sentence
        annot_idx = 0  # index of annotation that is within current sentence

        res = []
        M = len(self.annotations)
        for i, sentence in enumerate(all_sentences):
            # The sentence splitter isn't perfect. It recognizes "i.v." as a
            # sentence. Since there can be multiple instances of "sentences"
            # like "i.v." (e.g., PMID 10840460), we need to make sure that
            # we are checking for the first instance starting at the current
            # position (since find always finds the first instance otherwise).
            assert full_text.find(sentence, sent_idx) == sent_idx, (
                "PMID {0} sentence '{1}' does not match text!".format(self.pmid,
                                                                      sentence))

            sent_stop = sent_idx + len(sentence)

            start_annot = annot_idx
            while annot_idx < M and self.annotations[annot_idx].stop <= sent_stop:
                annot_idx += 1

            # should be one past
            res.append(Sentence(self.pmid, i, sentence,
                sent_idx, sent_stop, self.annotations[start_annot : annot_idx]))

            sent_idx += len(sentence) + 1  # all sentences separated by one space

        return res

    def _classify_relations(self):
        """
        Takes all the possible relations for this abstract and
        splits them into three mutually exclusive groups:
            1. CID relations, which are sentence bound
            2. Non-CID, sentence-bound relations
            3. Relations which are not sentence bound (head to abstract?) TODO
        """
        all_rels = self._get_all_possible_relations()

        cid_rels = set()
        sentence_non_cid_rels = set()
        for sentence in self.sentences:
            cid_rels |= sentence.possible_relations[True]
            sentence_non_cid_rels |= sentence.possible_relations[False]

        sentence_non_cid_rels -= cid_rels

        not_sent_bound_rels = all_rels - cid_rels - sentence_non_cid_rels

        assert cid_rels.is_disjoint(sentence_non_cid_rels)
        assert cid_rels.is_disjoint(not_sent_bound_rels)
        assert not_sent_bound_rels.is_disjoint(sentence_non_cid_rels)

        assert (len(self.chemicals) * len(self.diseases)
            == len(cid_rels | sentence_non_cid_rels | not_sent_bound_rels))

        possible_relations = {
            "CID": cid_rels,
            "sentence_non_CID": sentence_non_cid_rels,
            "not_sentence_bound": not_sent_bound_rels
        }
        return possible_relations

    def _has_relation(self, potential_relation):
        """
        Checks if the provided possible relation object matches any of the
        gold standard relations for this paper.

        Note:
            It is not possible to use a set to do the checking
            operation here, because it is not possible to make
            the hashes of two objects the same when they are
            defined by be equal by the overridden equals operator
            for Relation objects.

            This solution is slow, but at least it's correct.
        """
        return potential_relation in self.gold_relations


# TODO, add which sentence it is found in annotation? maybe

class Concept(models.Model):
    document = models.ForeignKey(Document)
    uid = models.TextField(blank=False)
    stype = models.TextField(blank=False)
    text = models.TextField(blank=False)
    start = models.IntegerField()
    stop = models.IntegerField()

    def __unicode__(self):
        return self.uid

class Sentence(models.Model):
    paper = models.ForeignKey(Paper)
    idx = models.TextField(blank=False)
    uid = models.TextField(blank=False)
    text = models.TextField(blank=False)
    start = models.IntegerField()
    stop = models.IntegerField()
    annotations = models.TextField(blank=False)

    def __unicode__(self):
        return self.text


class Relation(models.Model):
    document = models.ForeignKey(Document)
    relation = models.TextField(blank=False)
    concept1_id = models.ForeignKey(Concept, related_name='concept1')
    concept2_id = models.ForeignKey(Concept, related_name='concept2')

    automated_cid = models.BooleanField(default=False)
    # TODO if CID relation automatically determined, then apply
    # relation choice to auto
    if automated_cid == True:
        relation = "auto"

    def __unicode__(self):
        return self.relation


class Answer(models.Model):
    """
    Class where the user can provide their response to whether or not chemicals
    and diseases are related. TODO: if there is an auto CID, then don't have
    the user provide Answers to the relations!!
    """
    relation = models.ForeignKey(Relation)
    # TODO, this might get updated from Toby's code in the future
    relation_type = models.TextField(blank=True)
    # user confidence order 1 to 4 (where 1 is not confident and 4 is confident)
    # This value is used in scoring later

    username = models.TextField(blank=True)


"""
class Relation(models.Model):
    paper = models.ForeignKey(Paper) # Paper object is the foreign key
    pmid = models.TextField(blank=False)
    chemical_id = models.TextField(blank=False)
    disease_id = models.TextField(blank=False)

    def __unicode__(self):
        return self.chemical_id, self.disease_id
"""


"""

class Section(models.Model):
    SECTION_KIND_CHOICE = (
        ('t', 'title')
        ('a', 'abstract')
    )

    # kind of section it is and whether or not it has BOTH chemical and disease combinations
    # if not, then don't use it.
    kind = models.CharField(max_length=1, choices=SECTION_KIND_CHOICE)
    text = models.TextField(blank=True)

    # need the span of the sentences here for reference? numbers TODO (which program?)
    span_start = models.IntegerField()
    span_end = models.IntegerField()

    # at

    updated = models.DateTimeField(auto_now=True)
    created = models.DateTimeField(auto_now_add=True)

    # section relates back to DOCUMENT
    document = models.ForeignKey(Document)

    # TODO add more methods (look at examples)

"""
