from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Count

from mark2cure.document.models import Activity

from datetime import datetime
import datetime, random, re, nltk


def experiment_routing(user, n_count, gm_occurance = 4, k_max = 3):
    user_profile = user.userprofile

    # (TODO) Pretect from out of range & smarter gold injection
    gm_docs = [282, 310, 363, 326, 322]
    if n_count % gm_occurance == 0:
        gm_index = (n_count / gm_occurance) - 1
        if (gm_index+1) > len(gm_docs):
            return gm_docs[gm_index]

    if user_profile.mturk:
        prev_docs = Activity.objects.filter(user=user, experiment=settings.EXPERIMENT).values('document__pk', flat=True).all()
    else:
        prev_docs = Activity.objects.filter(user=user).values_list('document__pk', flat=True).all()

    '''
      I need to figure out which of the current experiment
      documents are still available for work to be done on
      1. Array of all document_id options
      2. Array of all documents already done by that user
      3. Array of all documents from options which have already been completed their max times
      4. Array of (1 - 2) - 3
    '''
    experiment_docs = [2787, 3195, 3357, 2726, 2637, 3030, 3203, 3314, 3077, 2369, 2394, 3003, 3567, 3166, 3177, 2152, 2661, 2236, 2193, 2878]
    for x in prev_docs:
      experiment_docs.remove(x)

    activities = Activity.objects.filter(document__pk__in=experiment_docs, task_type='cr', submission_type='gm', experiment=None).values('document').annotate(Count('document'))
    experiment_docs_completed = [item['document'] for item in activities if item['document__count'] >= k_max]

    for x in experiment_docs_completed:
        experiment_docs.remove(x)


    random.shuffle(experiment_docs)
    return experiment_docs[0]


def get_mturk_account(worker_id):
    u, created = User.objects.get_or_create(username=worker_id)
    if created:
        u.set_password('')
        profile = u.profile
        profile.email_notify = False
        profile.mturk = True
        profile.save()
        u.save()

    return u



def get_timezone_offset(tz):
    offset = datetime.now(tz).strftime('%z')

    sign = 1
    if offset[0] == '-':
        sign = -1

    return sign * float(offset[1:]) / 100

'''
  Text Utility Functions

'''

def is_negative(text):
    return False


def text_has_link(text):
    pattern = re.compile("(https?://|a href)")
    match = pattern.findall(text, re.I)
    return True if len(match) else False


def text_has_any_words(text, words):
    for word in text.split():
        if word in words:
            return True

    return False


def text_has_any_string(text, strings):
    for string in strings:
        if string in text:
            return True

    return False


def text_has_any_pattern(text, patterns):
    combined = "(" + "|".join(patterns) + ")"
    if re.match(combined, text):
        return True
    return False


'''
  Language Utilities
'''

def tokenize(text):
    try:
        return nltk.word_tokenize(text)
    except Exception, e:
        print e
        return text.split()

def pos_tag(text):
    tokens = tokenize(text)
    try:
        return nltk.pos_tag(tokens)
    except Exception:
        tags = []
        for word in tokens:
            tags.append((word, '???'))

        return tags


def get_proper_nouns(text):
    pos_tags = pos_tag(text)
    return [word for word,pos in pos_tags if pos == 'NNP']


def language_clean(text):
    text = text.strip()
    # text.delete!("\"") if text.count("\"").odd? # derp, we have an unmatched quote? -- let's remove ALL the quotes!
    # text.gsub!(/\A\)[\.\s]/, "") # dumbass parens
    # text = text.strip.strip_html.squeeze(" ").gsub(/\.{4,}/,
    # "...").gsub(/(\r|\n|\t)/," ").gsub(/\?{2,}/, "?").gsub(/\!{2,}/, "!")
    return text


def clean_spaces(text):
    return ' '.join(text.split()).strip()


def is_question(text):
    return '?' in text.strip()[-3:]


def is_respondable(text):
    return True


