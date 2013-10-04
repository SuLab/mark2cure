# -*- coding: utf-8 -*-
"""
    manage
    ~~~~~~

    Manager module
"""

from flask.ext.script import Manager

from mark2cure.api import create_app
from mark2cure.manage import *

manager = Manager(create_app())
manager.add_command('heatmap', Heatmap())
manager.add_command('annotate', Annotate())
manager.add_command('create', Create())
manager.add_command('gold', SolidGold())
manager.add_command('compare', Compare())
manager.add_command('turk', Turk())

if __name__ == "__main__":
  manager.run()
