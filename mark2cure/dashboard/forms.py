from django import forms

from ..common.models import Group

import re


class GroupForm(forms.ModelForm):
    pmids = forms.CharField(widget = forms.Textarea)

    def clean_pmids(self):
        pmids_text = self.cleaned_data['pmids']
        # PMIDs are all 8 ints
        pmids_arr = list(set(re.findall('\d{8}', pmids_text)))
        return pmids_arr

    class Meta:
        model = Group
        fields = ['name', 'stub', 'description',
                  'order', 'enabled', 'pmids']



