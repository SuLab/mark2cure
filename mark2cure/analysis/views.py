from django.template import RequestContext
from django.shortcuts import render_to_response
from django.shortcuts import render
from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.db.models import Count
from django.conf import settings


from mark2cure.document.models import *
from mark2cure.analysis.utils import clean

import json

@login_required
def experiment_details(request, exp_id):

    total_activities = Activity.objects.filter(task_type = 'cr', experiment = exp_id).exclude(user__userprofile__ignore = True).count()
    total_experimental_activities = Activity.objects.filter(task_type = 'cr', experiment = exp_id).exclude(submission_type = 'gm', user__userprofile__ignore = True).count()

    workers = Activity.objects.filter(
        task_type = 'cr',
        experiment = exp_id).exclude(user__userprofile__ignore = True).values('user', 'user__username').annotate(Count('user')).order_by('-user__count')

    # (TODO) multiple excludes?
    documents = Activity.objects.filter(
        task_type = 'cr',
        experiment = exp_id).exclude(submission_type = 'gm').values('document', 'document__document_id').annotate(Count('document')).order_by('-document__count')

    gold_documents = Activity.objects.filter(
        task_type = 'cr',
        submission_type = 'gm',
        experiment=exp_id).exclude(user__userprofile__ignore = True).values('document', 'document__document_id').annotate(Count('document')).order_by('-document__count')

    for doc in documents: doc['comment_count'] = Comment.objects.filter(document__pk = doc['document']).count()
    for doc in gold_documents: doc['comment_count'] = Comment.objects.filter(document__pk = doc['document']).count()

    f_scores = Activity.objects.filter(
        task_type = 'cr',
        submission_type = 'gm',
        experiment = exp_id).exclude(user__userprofile__ignore = True).values_list('f_score', flat = True).order_by('-f_score')

    f_scores_non_first = Activity.objects.filter(
        task_type = 'cr',
        submission_type = 'gm',
        experiment = exp_id).exclude(document_id__in = [869, 956, 1018, 520]).values_list('f_score', flat = True).order_by('-f_score')


    avg_f = reduce(lambda x, y: x + y, f_scores) / len(f_scores)


    return render_to_response('analysis/experiment_details.jade', {
      'exp_id' : exp_id,
      'total_activities' : total_activities,
      'total_experimental_activities': total_experimental_activities,

      'cost' : total_activities*.06,
      'avg_f' : avg_f,

      'workers' : workers,
      'documents' : documents,
      'gold_documents' : gold_documents,
      }, context_instance=RequestContext(request))




@login_required
def network(request):
    node_list = []
    link_list = []

    # (TODO) Refine to X # of articles, X days ago, ...
    # annotations = Annotation.objects.values('view__section__document', 'text').filter(view__user = request.user).distinct().all()
    annotations = Annotation.objects.filter(view__user = request.user).all()

    # Saving the unique arrays of annotations and documents
    ann_arr = list(set( [clean(ann.text)  for ann in annotations] ))
    doc_arr = list(set( [ann.view.section.document     for ann in annotations] ))

    # Put the uniq docs at the top of the list then uniq annotation terms after the documents
    [node_list.append({"name": doc.title, "group": 1})  for doc in doc_arr]
    [node_list.append({"name": ann, "group": 2})        for ann in ann_arr]

    # For each of the documents with annotations
    for doc_idx, doc in enumerate( doc_arr ):
      # Collect all of the cleaned annotations for that document by the current user
      doc_annotations = Annotation.objects.filter(view__section__document = doc, view__user = request.user).all()
      doc_anns = [clean(ann.text) for ann in doc_annotations]

      for ann in list(set(doc_anns)):
        # For this document, connect the doc to the annotation and weight appropriately
        link_list.append({ "source"  : doc_idx,
                           "target"  : len(doc_arr) + ann_arr.index(ann),
                           "value"   : doc_anns.count(ann) })

        return HttpResponse(json.dumps({'nodes' : node_list, 'links' : link_list}), content_type="application/json")

