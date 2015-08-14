from django.template.response import TemplateResponse


def disease_marking(request):
    return TemplateResponse(request, 'instructions/disease-marking.jade')


def gene_marking(request):
    return TemplateResponse(request, 'instructions/gene-marking.jade')


def treatment_marking(request):
    return TemplateResponse(request, 'instructions/treatment-marking.jade')


def relation_identification(request):
    return TemplateResponse(request, 'instructions/relation-identification.jade')
