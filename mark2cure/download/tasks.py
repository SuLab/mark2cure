from ..document.models import Document


def group_export(request, document_pks, er=True, rel=True):
    writer = Document.objects.as_writer(documents=document_pks)

    '''
    if er:
        for doc in group.get_documents():
            doc_writer = doc.as_writer()
            if selection_type == 'user':
                doc_df = doc.as_er_df_with_user_annotations()
            else:
                doc_df = doc.as_er_df_with_pubtator_annotations()
            # Protection isn't needed b/c this is the raw output for analysis.
            doc_df = clean_df(doc_df, overlap_protection=False)
            # Convert DF table into BioC Document
            doc_writer = apply_er_annotations(doc_writer, doc_df)
            writer.collection.add_document(doc_writer.collection.documents[0])
    '''

    if rel:
        doc_df = Document.objects.relation_df(documents=document_pks)

    # doc_writer = apply_rel_annotations(doc_writer, doc_df)
    # writer.collection.add_document(doc_writer.collection.documents[0])

