# -*- coding: utf-8 -*-
"""
    tests.api.text
    ~~~~~~~~~~~~~~~~~~~~

"""

from . import Mark2CureApiTestCase


class TextApiTestCase(Mark2CureApiTestCase):

    def test_get_current_user(self):
        self.assertOkJson('{}')
