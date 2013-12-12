from mark2cure.common.utils import Turk

turk = Turk()

# turk.disable_all()
# turk.make_qualification_test()
turk.hit_for_document(112)



# documents = db.session.query(Document).\
#     filter_by(source = 'NCBI_corpus_development').\
#     all()
# doc_ids = [doc.id for doc in documents]
# random.shuffle(doc_ids)
# for doc_id in doc_ids:
#   self.hit_for_document( doc_id )
