from mark2cure.document.models import *

EXPERIMENT_DOCS = [363, 290, 374, 296, 368, 311, 349, 284, 322, 355, 307, 371, 280, 314, 334, 364, 319, 305, 341, 303, 310, 293, 372, 333, 320, 373, 291, 302, 326, 351, 325, 367, 327, 294, 369, 289, 336, 318, 283, 313, 350, 304, 362, 361, 357, 352, 282, 299, 348, 344]

print "Non GM ID", "Original GM ID"

for doc_id in EXPERIMENT_DOCS:
    old_gm_doc = Document.objects.filter(pk = doc_id).first()
    non_gm_doc = Document(document_id = old_gm_doc.document_id, title = old_gm_doc.title, source = 'NCBI_corpus_development-NON-GM')
    non_gm_doc.save()
    old_gm_sec = old_gm_doc.section_set.filter(kind='t').first()
    non_gm_sec = Section(kind = 't', document = non_gm_doc, text = old_gm_sec.text)
    non_gm_sec.save()
    old_gm_sec = old_gm_doc.section_set.filter(kind='a').first()
    non_gm_sec = Section(kind = 'a', document = non_gm_doc, text = old_gm_sec.text)
    non_gm_sec.save()
    print non_gm_doc.pk, old_gm_doc.pk

