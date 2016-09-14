from django.template.response import TemplateResponse
from ..common.models import Group
from ..task.relation.models import RelationGroup
from django.views.decorators.http import require_http_methods
from django.shortcuts import redirect

from .models import Download


def home(request):
    groups = [(x[0], x[1], 'er') for x in Group.objects.values_list('pk', 'name')] + [(x[0], x[1], 'rel') for x in RelationGroup.objects.values_list('pk', 'name')]
    ctx = {
        'groups': groups,
        'downloads': Download.objects.all()
    }
    return TemplateResponse(request, 'download/home.jade', ctx)


@require_http_methods(['POST'])
def start_export(request):
    task_type = request.POST.get('task_type')
    group_pk = request.POST.get('group_pk')

    if task_type == 'er':
        group = Group.objects.get(pk=group_pk)
        docs = group.get_documents()
    else:
        group = RelationGroup.objects.get(pk=group_pk)
        docs = group.documents.all()

    d = Download.objects.create(
        task_er=task_type == 'er',
        task_rel=task_type == 'rel',
    )
    d.documents = docs.values_list('pk', flat=True)
    d.save()

    return redirect('download:home')

