# Tong Shu Li
# First written: 2015-07-02
# Last updated: 2015-07-30
"""
Data models for BioCreative V task 3.

In order to make processing the BioCreative
data easier, these data models were created
in order to simplify the tasks of:
1. Parsing the gold standard
2. Concept highlighting
3. Concept co-occurrence detection
4. Relationship verification
5. Work unit generation
"""
from collections import defaultdict

from lingpipe.file_util import read_file
from lingpipe.split_sentences import split_abstract


def is_MeSH_id(uid):
    return len(uid) == 7 and uid[0] in ["C", "D"]


class Annotation:
    """
    A single mention of a concept in a piece of text.
    Annotation positions are indexed to the abstract.
    """
    def __init__(self, uid, stype, text, start, stop):
        if ":" in uid:
            self.uid_type, self.uid = uid.split(':')
        else:
            self.uid_type = "unknown"
            self.uid = uid

        self.stype = stype.lower()
        assert self.stype in ["chemical", "disease"]

        self.text = text
        self.start = int(start)
        self.stop = int(stop)
        assert self.start < self.stop, "Annotation {0} indicies reversed!".format(self.uid)
        assert len(text) == self.stop - self.start, "Annotation {0} length mismatch!".format(self.uid)

    def __repr__(self):
        return "<{0}>: '{1}'({2}) {3}-{4}".format(self.__class__.__name__,
            self.text, self.stype, self.start, self.stop)

    def __cmp__(self, other):
        """
        When sorting Annotation objects, make sure they
        are in ascending order of starting position.
        """
        if hasattr(other, "start"):
            return self.start.__cmp__(other.start)


class Sentence:
    """
    A single sentence from a paper.

    Sentence-bound relationships are generated
    by this object.

    The set of chemical-disease relations that need to
    be made into a work unit is the non-CID relations
    minus the CID relations true at the abstract level.
    """
    def __init__(self, pmid, idx, text, start, stop, annotations):
        self.pmid = int(pmid)
        self.uid = "{0}_{1}".format(pmid, idx)
        self.text = text
        self.start = int(start)
        self.stop = int(stop)
        assert self.start < self.stop, "Sentence {0} indicies reversed!".format(self.uid)
        assert len(text) == self.stop - self.start, "Sentence {0} length mismatch!".format(self.uid)

        # a list of the concept annotations within this sentence
        self.annotations = annotations

        # generate the list of CID and non-CID relations bound to this sentence
        self.poss_relations = self.classify_relations()

    def __repr__(self):
        return "<{0}>: PMID: {1} '{2}'({3}-{4})\nAnnotations: {5}\n".format(
            self.__class__.__name__, self.pmid, self.text,
            self.start, self.stop, self.annotations)

    def is_CID_relation(self, chemical, disease):
        """
        Given two Annotations within this sentence,
        determine if the pair of annotations follows the
        CID structure.
        """
        return ((chemical.stop < disease.start)
            and (disease.start - chemical.stop <= 15)
            and ("induce" in self.text[chemical.stop - self.start:
                disease.start - self.start].lower())
        )

    def classify_relations(self):
        """
        This function generates all the unique chemical-disease
        identifier pairs of annotations contained within this
        sentence, and classifies them into two groups:
            1. Those which follow a '[chemical]-induced [disease]' (CID)
                structure.
            2. Those which do not follow the CID structure.

        The two sets of identifier pairs are mutually exclusive.
        """
        all_relations = defaultdict(set)
        for annot_A in self.annotations:
            if annot_A.stype == "chemical" and annot_A.uid != "-1":
                for annot_B in self.annotations:
                    if annot_B.stype == "disease" and annot_B.uid != "-1":
                        rel_type = self.is_CID_relation(annot_A, annot_B)
                        all_relations[rel_type].add((annot_A.uid, annot_B.uid))

        """
        In cases where we have a sentence with the following annotations:
        D C D (where D = disease and C = chemical), we see that the first
        instance of C and D is not in the CID format, and will get added
        to the non-CID set. However, the second instance is in the CID
        format, and gets added to the CID set. To make sure the two sets
        are mutually exclusive, we need to subtract the CID set from the
        larger non-CID set.
        """
        all_relations[False] -= all_relations[True]
        assert all_relations[False].isdisjoint(all_relations[True])
        return all_relations


class Relation:
    """
    A single chemical-induced disease relationship.
    Used to compare identifier pairs against the gold.
    """
    def __init__(self, pmid, chemical_id, disease_id):
        self.pmid = int(pmid)
        self.uid = "{0}:{1}-{2}".format(pmid, chemical_id, disease_id)

        assert chemical_id != "-1" and disease_id != "-1", "Relation {0} has bad ids.".format(self.uid)
        assert is_MeSH_id(chemical_id), "Relation {0} has bad chemical id.".format(self.uid)

        # disease ids can be complex ones joined together
        # represent the disease ids as a set
        disease_id = set(disease_id.split('|'))

        for uid in disease_id:
            assert is_MeSH_id(uid), "Relation {0} has bad disease id.".format(self.uid)

        self.chemical_id = chemical_id
        self.disease_id = disease_id

    def __eq__(self, other):
        """
        Equal if the chemical ids match exactly and
        at least one of the disease ids are shared
        between the two relations.

        This is because the gold relations only use
        a pair of single MeSH ids, despite the fact
        that the annotations use complexed MeSH ids.

        WARNING:
            Defining the equals function in this manner
            breaks transitivity. That is, if we have
            three objects A, B, and C, then if
                A == B and B == C, then
                A == C IS NOT TRUE!!
        """
        if isinstance(other, self.__class__):
            return (self.chemical_id == other.chemical_id
                and len(self.disease_id & other.disease_id) > 0)

        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return "<{0}>: {1}".format(
            self.__class__.__name__, self.uid)


class Paper:
    """
    A single academic publication.
    Contains:
        1. The PubMed identifier.
        2. The title as a string.
        3. The abstract as a string.
        4. A list of all chemical and disease annotations
            in the title and abstract sorted in increasing
            order of starting index.
        5. A potentially empty list of gold standard CID
            relations.
        6. A set of unique chemical identifiers.
        7. A set of unique disease identifiers.
        8. A list of Sentences containing both the title
            and body of the abstract. The first sentence
            is the title. Each sentence contains the
            annotations and relations constrained to that
            particular sentence.
        9. A set of all the potential chemical-disease
            relations grouped into three mutually exclusive
            categories:
            - CID relations
            - Non-CID sentence-bound relations
            - Non-sentence bound relations
            The sum of relations in all three groups should
            equal the number of unique chemical IDs times
            the number of unique disease IDs.
    """
    def __init__(self, pmid, title, abstract, annotations, gold_relations=[]):
        self.pmid = int(pmid)
        self.title = title
        self.abstract = abstract

        self.annotations = sorted(annotations)
        assert self.has_correct_annotations()

        self.gold_relations = gold_relations  # may be empty when not parsing gold

        self.chemicals, self.diseases = self.get_unique_concepts()

        # split sentences and generate sentence-bound relations
        self.sentences = self.split_sentences()
        self.poss_relations = self.classify_relations()

    def __repr__(self):
        return ("<{0}>: PMID {1}. {2} annotations, {3} gold relations\n"
            "{4} unique chemical ids, {5} unique disease ids\n"
            "{6} sentences".format(self.__class__.__name__,
            self.pmid, len(self.annotations), len(self.gold_relations),
            len(self.chemicals), len(self.diseases), len(self.sentences))
        )

    def has_correct_annotations(self):
        """
        Checks that the paper's annotations match the
        stated positions in the text.
        """
        text = "{0} {1}".format(self.title, self.abstract)
        for annotation in self.annotations:
            assert text[annotation.start: annotation.stop] == annotation.text, (
                "Annotation {0} in PMID {1} does not match the text.".format(annotation, self.pmid))

        return True

    def get_unique_concepts(self):
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

    def get_all_possible_relations(self):
        """
        Returns all possible unique drug-disease combinations
        as a set for this paper.
        """
        return {(chemical_id, disease_id) for chemical_id in self.chemicals for disease_id in self.diseases}

    def split_sentences(self):
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
                "PMID {0} sentence '{1}' does not match text!".format(self.pmid, sentence))

            sent_stop = sent_idx + len(sentence)

            start_annot = annot_idx
            while annot_idx < M and self.annotations[annot_idx].stop <= sent_stop:
                annot_idx += 1

            # should be one past
            res.append(Sentence(self.pmid, i, sentence,
                sent_idx, sent_stop, self.annotations[start_annot: annot_idx]))

            sent_idx += len(sentence) + 1  # all sentences separated by one space

        return res

    def classify_relations(self):
        """
        Takes all the possible relations for this abstract and
        splits them into three mutually exclusive groups:
            1. CID relations, which are sentence bound
            2. Non-CID, sentence-bound relations
            3. Relations which are not sentence bound
        """
        all_rels = self.get_all_possible_relations()

        cid_rels = set()
        sentence_non_cid_rels = set()
        for sentence in self.sentences:
            cid_rels |= sentence.poss_relations[True]
            sentence_non_cid_rels |= sentence.poss_relations[False]

        sentence_non_cid_rels -= cid_rels

        not_sent_bound_rels = all_rels - cid_rels - sentence_non_cid_rels

        assert cid_rels.isdisjoint(sentence_non_cid_rels)
        assert cid_rels.isdisjoint(not_sent_bound_rels)
        assert not_sent_bound_rels.isdisjoint(sentence_non_cid_rels)

        assert (len(self.chemicals) * len(self.diseases)
            == len(cid_rels | sentence_non_cid_rels | not_sent_bound_rels))

        poss_relations = {
            "CID": cid_rels,
            "sentence_non_CID": sentence_non_cid_rels,
            "not_sentence_bound": not_sent_bound_rels
        }
        return poss_relations

    def has_relation(self, potential_relation):
        """
        Checks if the provided possible Relationship object
        matches any of the gold standard relationships for this
        paper.

        Note:
            It is not possible to use a set to do the checking
            operation here, because it is not possible to make
            the hashes of two objects the same when they are
            defined by be equal by the overridden equals operator
            for Relation objects.

            This solution is slow, but at least it's correct.
        """
        return potential_relation in self.gold_relations


def parse_input(loc, fname, is_gold=True, return_format="list"):
    """
    Reads a given file and returns a list of Paper
    objects.
    """
    assert return_format in ["list", "dict"]
    if return_format == "list":
        papers = []
    else:
        papers = dict()

    counter = 0
    annotations = []
    relations = []
    for i, line in enumerate(read_file(fname, loc)):
        if len(line) == 0:
            # time to create the paper object
            if return_format == "list":
                papers.append(Paper(pmid, title, abstract, annotations, relations))
            else:
                papers[pmid] = Paper(pmid, title, abstract, annotations, relations)

            counter = 0
            annotations = []
            relations = []
        else:
            if 0 <= counter <= 1:
                vals = line.split('|')
                assert len(vals) == 3, "Title or abstract on line {i} is messed up!".format(i + 1)
            else:
                vals = line.split('\t')

            if counter == 0:
                assert vals[1] == "t", i + 1
                pmid = int(vals[0])
                title = vals[2]
            elif counter == 1:
                assert vals[1] == "a"
                assert int(vals[0]) == pmid
                abstract = vals[2]
            elif is_gold and len(vals) == 4:
                assert vals[1] == "CID"
                relations.append(Relation(pmid, vals[2], vals[3]))
            else:
                assert 6 <= len(vals) <= 7
                annotations.append(Annotation(vals[5], vals[4], vals[3], vals[1], vals[2]))

            counter += 1

    return papers
