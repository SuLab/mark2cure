#!/usr/bin/env python
import os
import sys

if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "mark2cure.settings")
    os.environ.setdefault("DJANGO_CONFIGURATION", "Development")

    from configurations import importer
    importer.install()

    from configurations.management import execute_from_command_line
    execute_from_command_line(sys.argv)
