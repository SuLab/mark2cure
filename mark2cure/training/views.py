from django.template.response import TemplateResponse
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.shortcuts import redirect

from ..api.views import get_training_dict


@login_required
def route(request):

    arr_of_tasks_completed_training = [task['task'] for task in get_training_dict(request.user.pk) if len(task['levels']) > 0 if all([lvl.get('completions') > 0 for lvl in task['levels']])]

    if len(arr_of_tasks_completed_training) == 2:
        return redirect('common:dashboard')
    else:
        return redirect(reverse('training:re'))

    # (TODO) Is there ever a time we'd use this [requires this training to be completed first!]
    # return redirect(reverse('training:ner'))


@login_required
def ner_home(request):
    '''Starting HTML for all NER Training Javascript
        (TODO) Not yet built
    '''
    return TemplateResponse(request, 'training/entity-recognition/home.html')


@login_required
def re_home(request):
    '''Starting HTML for all RE Training Javascript
    '''
    return TemplateResponse(request, 'training/relation/home.html')

