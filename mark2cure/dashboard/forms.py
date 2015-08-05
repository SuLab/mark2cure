from django import forms

from ..common.models import Group

import re


class GroupForm(forms.ModelForm):
    pmids = forms.CharField(widget = forms.Textarea)

    def clean_pmids(self):
        pmids_text = self.cleaned_data['pmids']
        # Any digit, no length checks (PMIDs are all over)
        pmids_arr = list(set(re.findall('\d+', pmids_text)))
        return pmids_arr

    class Meta:
        model = Group
        fields = ['name', 'stub', 'description',
                  'order', 'enabled', 'pmids']


