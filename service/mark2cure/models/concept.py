from . import db

class Concept(db.Model):
    id          = db.Column(db.Integer, primary_key=True)
    concept_id  = db.Column(db.Text)

    annotations = db.relationship('Annotation',   backref=db.backref('concept',  lazy='select'))

    def __init__(self, concept_id):
        self.concept_id = concept_id
