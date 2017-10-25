from django.template.response import TemplateResponse
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.shortcuts import redirect

from ..task.models import Task
from ..task.models import Level


@login_required
def route(request):
    # Prioritize relation training over Entity Recognition for routing
    if Level.objects.filter(user=request.user, task_type='re').exists():
        return redirect(reverse('training:re'))

    user_level = Level.objects.filter(user=request.user, task_type='ner').first().level
    if user_level <= 3:
        task = Task.objects.get(kind=Task.TRAINING, provides_qualification='4')
    else:
        task = Task.objects.filter(kind=Task.TRAINING, requires_qualification=user_level).first()

    if task:
        return redirect(task.meta_url)
    else:
        return redirect('common:dashboard')


def ner_home(request):
    '''Starting HTML for all NER Training Javascript
        (TODO) Not yet built
    '''
    return TemplateResponse(request, 'training/entity-recognition/home.html')


def re_home(request):
    '''Starting HTML for all RE Training Javascript
    '''
    return TemplateResponse(request, 'training/relation/home.html')

