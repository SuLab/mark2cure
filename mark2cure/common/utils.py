from django.conf import settings
from datetime import datetime

import random, re, nltk

def get_timezone_offset(tz):
    offset = datetime.now(tz).strftime('%z')

    sign = 1
    if offset[0] == '-':
        sign = -1

    return sign * float(offset[1:]) / 100

'''
  Text Utility Functions

  this whole thing is fucked

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
