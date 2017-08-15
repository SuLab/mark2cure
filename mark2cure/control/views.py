from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import user_passes_test
from django.contrib.auth.decorators import login_required
from django.template.response import TemplateResponse
from django.shortcuts import get_object_or_404, redirect
from django.core.urlresolvers import reverse
from django.contrib.auth.models import User
from django.http import HttpResponse, HttpResponseRedirect

from .forms import GroupForm
from ..analysis.models import Report
from ..userprofile.models import UserProfile
from ..document.models import Document, Pubtator
from ..document.tasks import get_pubmed_document

from ..common.models import Group
from ..task.models import Task, UserQuestRelationship

import pandas as pd
import itertools
import datetime
import uuid


def dataframe_view(request, df, format_type, ctx={}, template='control/dataframe_base.jade'):
    if format_type == "download":
        return HttpResponse(df.to_csv(), content_type='text/csv')

    if format_type == "csv":
        return HttpResponse(df.to_csv())

    elif format_type == "json":
        return HttpResponse(df.to_json(), content_type='application/json')

    elif format_type == "string":
        return HttpResponse(df.to_string)

    else:
        stnd_ctx = {
            'dataframe': df,
            'dataframe_html': df.to_html(classes='table table-striped table-condensed'),
            'view_name': request.resolver_match.namespace + ':' + request.resolver_match.url_name
        }
        stnd_ctx.update(ctx)
        return TemplateResponse(request, template, stnd_ctx)


@login_required
@user_passes_test(lambda u: u.is_staff)
def user_training(request, format_type="html"):
    training = Task.objects.filter(kind=Task.TRAINING).all()

    arr = []
    for user in User.objects.all():
        row = [user.username, user.pk]
        [row.append(UserQuestRelationship.objects.filter(task=t, user=user, completed=True).exists()) for t in training]
        arr.append(row)

    columns = itertools.chain(["Username", "User PK"], training.values_list('pk', flat=True))
    df = pd.DataFrame(arr, columns=list(columns))
    return dataframe_view(request, df, format_type)


@login_required
@user_passes_test(lambda u: u.is_staff)
def user_quest_availability(request, format_type="html"):
    # (TODO) Broken with task changes and notions of quest completions with Relation app
    return

    arr = []
    users = User.objects.exclude(userprofile__isnull=True).all()
    for u in users.values('pk', 'username'):
        try:
            profile = UserProfile.objects.get(user__pk=u['pk'])
        except Exception as e:
            pass

    df = pd.DataFrame(arr, columns=['Username', 'Quest Count', 'Available Quests'])
    return dataframe_view(request, df, format_type)


@login_required
@user_passes_test(lambda u: u.is_staff)
@require_http_methods(['GET', 'POST'])
def group_create(request):

    if request.method == 'POST':
        group_form = GroupForm(instance=None, data=request.POST or None)
        group_uuid = str(uuid.uuid1())

        if group_form.is_valid():
            group = group_form.save()
            pmids = group_form.cleaned_data.get('pmids', [])
            get_pubmed_document.delay(pmids, source=group_uuid, group_pk=group.pk)
            return redirect(reverse('control:group_list'))

    if request.method == 'GET':
        group_form = GroupForm(instance=None, data=request.POST or None)
        ctx = {
            'form': group_form
        }
        return TemplateResponse(request, 'control/group_create.jade', ctx)


@login_required
@user_passes_test(lambda u: u.is_staff)
def group_list(request):
    ctx = {
        'groups': Group.objects.all().order_by('-order'),
    }
    return TemplateResponse(request, 'control/group_list.jade', ctx)


@login_required
@user_passes_test(lambda u: u.is_staff)
def group_report(request, pk, format_type='html'):
    report = get_object_or_404(Report, pk=pk)
    ctx = {
        'model_pk': report.pk
    }
    return dataframe_view(request, report.dataframe, format_type, ctx, template='control/dataframe_pk_base.jade')


@login_required
@user_passes_test(lambda u: u.is_staff)
def group_read(request, pk):
    group = get_object_or_404(Group, pk=pk)
    document_quest_relationships = group.total_documents().select_related('task', 'document')

    ctx = {
        'group': group,
        'document_quest_relationships': document_quest_relationships
    }
    return TemplateResponse(request, 'control/group.jade', ctx)


@login_required
@require_http_methods(['POST'])
@user_passes_test(lambda u: u.is_staff)
def pubtator_actions(request, pk):
    pubtator = get_object_or_404(Pubtator, pk=pk)
    doc = pubtator.document

    if request.method == 'POST':
        pubtator.delete()
        doc.run_pubtator()

    return redirect(reverse('control:document', kwargs={'pk': doc.pk}))


@login_required
@require_http_methods(['POST'])
@user_passes_test(lambda u: u.is_staff)
def document_pubtator_actions(request, pk):
    doc = get_object_or_404(Document, pk=pk)

    if request.method == 'POST':
        doc.pubtators.all().delete()
        doc.run_pubtator()

    return HttpResponseRedirect(request.META.get('HTTP_REFERER'))


@login_required
@user_passes_test(lambda u: u.is_staff)
def document_read(request, pk):
    doc = get_object_or_404(Document, pk=pk)
    ctx = {
        'doc': doc
    }
    return TemplateResponse(request, 'control/doc.jade', ctx)


@login_required
@user_passes_test(lambda u: u.is_staff)
def home(request):
    today = datetime.datetime.now()
    users_online = {
        'hour': UserProfile.objects.filter(last_seen__gte=today - datetime.timedelta(hours=1)).count(),
        'day': UserProfile.objects.filter(last_seen__gte=today - datetime.timedelta(days=1)).count(),
        'week': UserProfile.objects.filter(last_seen__gte=today - datetime.timedelta(days=7)).count(),
        'month': UserProfile.objects.filter(last_seen__gte=today - datetime.timedelta(days=30)).count()
    }
    ctx = {
        'users_online': users_online
    }
    return TemplateResponse(request, 'control/home.jade', ctx)

