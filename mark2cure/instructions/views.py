from django.template.response import TemplateResponse


def disease_marking(request):
    return TemplateResponse(request, 'instructions/disease-marking.jade')


def gene_marking(request):
    return TemplateResponse(request, 'instructions/gene-marking.jade')


def treatment_marking(request):
    return TemplateResponse(request, 'instructions/treatment-marking.jade')


def gene_disease_relation(request):
    return TemplateResponse(request, 'instructions/gene-disease-relation.html')


def drug_disease_relation(request):
    return TemplateResponse(request, 'instructions/drug-disease-relation.html')


def gene_drug_relation(request):
    return TemplateResponse(request, 'instructions/gene-drug-relation.html')


# (TODO) this will be for a later date when we add more relations.
def gene_gene_relation(request):
    pass


def drug_drug_relation(request):
    pass


def disease_disease_relation(request):
    pass
