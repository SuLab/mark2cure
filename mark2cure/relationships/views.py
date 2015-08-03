import datetime
import re

# learn shortcuts
from django.shortcuts import get_object_or_404, render
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.views import generic
from django.utils import timezone


from .models import Paper
from .tasks import parse_input

import os
# TODO check imports
### add here if not completed TODO (remove from list if user already answered)
"""
Toby's parse_input creates new Paper objects using all of the
inputs to the Paper class.

"""


#def parse_input(location, fname, is_gold = True, return_format = "list"):


class IndexView(generic.ListView):
    template_name = 'relationships/index.html'

    # call the paper parsing function in tasks.py
    #TODO
    # make new variable for the "view" to use in a print statement TODO
    #papers = parse_input(os.getcwd(),"CDR_small.txt")

    #toby_test = "HELLO"

    def get_queryset(self):
        return Paper.objects.all()

class DetailView(generic.DetailView):
    model = Paper
    template_name = 'relationships/detail.html'

    def get_queryset(self):
        return Paper.objects.all()

class ResultsView(generic.DetailView):
    model = Paper
    template_name = 'relationships/results.html'

def vote(request, question_id):
    p = get_object_or_404(Paper, pk=Paper.pmid)
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
