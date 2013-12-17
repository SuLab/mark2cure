from mark2cure.common.utils import Turk
from mark2cure.document.models import Document
import random

turk = Turk()

turk.disable_all()
# turk.make_qualification_test()


################################
# Post the experiment documents
################################
documents = Document.objects.filter(source = 'NCBI_corpus_development').all()[:10]
doc_ids = [doc.id for doc in documents]
random.shuffle(doc_ids)
print doc_ids
for doc_id in doc_ids:
  turk.hit_for_document(doc_id)


