from django.template.response import TemplateResponse
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.shortcuts import redirect

from ..task.models import Task
from ..task.models import Level


# task = Task.objects.filter(kind=Task.TRAINING, provides_qualification=qualification_level).first()
# UserQuestRelationship.objects.create(task=task, user=user, completed=True)
# # Assign points to a the specific training level
# Point.objects.get_or_create(user=user,
#                             amount=task.points,
#                             content_type=ContentType.objects.get_for_model(task),
#                             object_id=task.id)
#
# Level.objects.create(user=user, task_type='ner', level=task.provides_qualification, created=timezone.now())


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
    return TemplateResponse(request, 'training/entity-recognition/home.html')


def re_home(request):
    return TemplateResponse(request, 'training/relation/home.html')

