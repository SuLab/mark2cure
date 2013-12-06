from mark2cure.document.models import Document

# # Get the gold bot account to make the db entries
# user = db.session.query(User).get(2)
# with open('../assets/NCBI_corpus/'+ path +'.txt','r') as f:
#   reader = csv.reader(f, delimiter='\t')
#   for num, title, text in reader:
#     text = text.replace('<category=', '<category type=')
#     plain_text = str(strip_tags(text, ['category']))
#
#     title = title.replace('<category=', '<category type=')
#     plain_title = str(strip_tags(title, ['category']))
#
#     soup = BeautifulSoup(text)
#
#     # Add the document to the database
#     doc = Document( int(num),
#                     plain_text,
#                     plain_title,
#                     datetime.datetime.utcnow(),
#                     path )
#     db.session.add(doc)
#     db.session.commit()
#
#     cats = set([cat.text for cat in soup.findAll('category')])
#     for cat in cats:
#       for m in re.finditer( cat, plain_text ):
#         # print( cat.text, cat['type'], m.start(), m.end() )
#         ann = Annotation( 0,
#                           'disease',
#                           cat,
#                           m.start(),
#                           len(cat),
#                           m.end(),
#                           user,
#                           doc,
#                           'gold_1.0',
#                           '',
#                           None );
#         db.session.add(ann)
#         # Save every document instead of once incase some doc crashes
#         db.session.commit()
