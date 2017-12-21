from django.template.response import TemplateResponse
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.shortcuts import redirect

from ..api.views import get_training_dict


@login_required
def route(request):

    res = get_training_dict(request.user.pk)
    print(res)
    # return redirect('common:dashboard')
    # return redirect(reverse('training:re'))
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

