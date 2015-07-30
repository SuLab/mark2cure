import datetime
import re

from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.views import generic
from django.utils import timezone


from .models import Choice, Question
###jennifer

### add here if not completed TODO (remove from list if user already answered)
def _create_question(question_text, days):
    time = timezone.now() + datetime.timedelta(days=days)
    return Question.objects.create(question_text=question_text,
                                   pub_date=time)

class IndexView(generic.ListView):
    template_name = 'relationships/index.html'
    context_object_name = 'latest_question_list'

    ### Small drug list, with less combinations:
    drug_list = ['ibuprofen','caffeine','slimfast', 'snickers bars']
    disease_list = ['obesity','headaches','fatigue','hypertension']

    for i in drug_list:
        for j in disease_list:
            new_question_text = "The sentence states that %s and %s:"%(i,j)
            if new_question_text not in str(Question.objects.all()):
                p = _create_question(new_question_text,0)   ###TODO remove publication date... this will become
                p.choice_set.create(choice_text="Are definitely associated (positive)")
                p.choice_set.create(choice_text="Are speculatively associated (speculative)")
                p.choice_set.create(choice_text="Are not associated (negative)")
                p.choice_set.create(choice_text="No claim of association made (false)")

    def get_queryset(self):

        """
        Return the last 100 published questions (not including those set to be
        published in the future).
        """
        return Question.objects.filter(pub_date__lte=timezone.now()).order_by('-pub_date')[:100]

class DetailView(generic.DetailView):
    model = Question
    template_name = 'relationships/detail.html'

    def get_queryset(self):
        """
        Excludes any questions that aren't published yet.
        """
        return Question.objects.filter(pub_date__lte=timezone.now())

class ResultsView(generic.DetailView):
    model = Question
    template_name = 'relationships/results.html'

def vote(request, question_id):
    p = get_object_or_404(Question, pk=question_id)
    try:
        selected_choice = p.choice_set.get(pk=request.POST['choice'])
    except (KeyError, Choice.DoesNotExist):
        # Redisplay the question voting form.
        return render(request, 'relationships/detail.html', {
            'question': p,
            'error_message': "Error: You didn't select one of the choices :(",
        })
    else:
        selected_choice.votes += 1
        selected_choice.save()
        def _create_question(question_text, days):
            time = timezone.now() + datetime.timedelta(days=days)
            return Question.objects.create(question_text=question_text,
                                           pub_date=time)
        if "positive" in str(selected_choice) or "speculative" in str(selected_choice): # Populates a new question when positive or speculative
            test_string_relationship = p.question_text[24:-1]
            r = _create_question("The relationship between"+test_string_relationship+" can be further described as...",0)
            replacement_text1 = re.sub("and", "may cause", test_string_relationship)
            replacement_text2 = re.sub("and", "is used in the treatment of", test_string_relationship)
            r.choice_set.create(choice_text=replacement_text1)
            r.choice_set.create(choice_text=replacement_text2)
            r.choice_set.create(choice_text="No extra information is given.")
            r.choice_set.create(choice_text="There is extra information that is not given as a choice here.")

        # Always return an HttpResponseRedirect after successfully dealing
        # with POST data. This prevents data from being posted twice if a
        # user hits the Back button.
        return HttpResponseRedirect(reverse('relationships:results', args=(p.id,)))
