import re

def clean(text):
  word = re.sub(r'\W+', ' ', text).lower().strip()
  return word[:-1] if word in ['cancers', 'chordomas'] else word


