from django.template.response import TemplateResponse
from ..common.models import Group
from ..task.relation.models import RelationGroup
from django.views.decorators.http import require_http_methods
from django.shortcuts import redirect

from .models import Download
from .tasks import group_export


def home(request):
    groups = [(x[0], x[1], 'er') for x in Group.objects.values_list('pk', 'name')] + [(x[0], x[1], 'rel') for x in RelationGroup.objects.values_list('pk', 'name')]
    ctx = {
        'groups': groups,
        'downloads': Download.objects.all()
    }
    return TemplateResponse(request, 'download/home.html', ctx)


@require_http_methods(['POST'])
def start_export(request):
    task_type = request.POST.get('task_type')
    group_pk = request.POST.get('group_pk')
    document_pks = request.POST.get('document_pks')
    ALL = 0
    NER = 1
    RE = 2

    if task_type == 'ner':
        group = Group.objects.get(pk=group_pk)
        docs = group.get_documents()
        group_export.apply_async(
            args=[list(docs.values_list('pk', flat=True))],
            kwargs={'export_type': NER},
            queue='mark2cure_downloads')

    elif task_type == 're':
        group = RelationGroup.objects.get(pk=group_pk)
        docs = group.documents.all()
        group_export.apply_async(
            args=[list(docs.values_list('pk', flat=True))],
            kwargs={'export_type': RE},
            queue='mark2cure_downloads')

    else:
        group_export.apply_async(
            args=[document_pks],
            kwargs={'export_type': ALL},
            queue='mark2cure_downloads')

    return redirect('download:home')

