from django.template.response import TemplateResponse


def disease_marking(request):
    return TemplateResponse(request, 'instructions/disease-marking.jade')


def gene_marking(request):
    return TemplateResponse(request, 'instructions/gene-marking.jade')


def treatment_marking(request):
    return TemplateResponse(request, 'instructions/treatment-marking.jade')


def gene_disease_relation(request):
    return TemplateResponse(request, 'instructions/gene-disease-relation.html')

def chemical_disease_relation(request):
    return TemplateResponse(request, 'instructions/chemical-disease-relation.html')

def gene_chemical_relation(request):
    # TODO
    pass

""" TODO this will be for a later date when we add more relations."""
def gene_gene_relation(request):
    pass

def chemical_chemical_relation(request):
    pass

def disease_disease_relation(request):
    pass
