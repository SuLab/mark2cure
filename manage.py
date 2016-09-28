#!/usr/bin/env python
import os
import sys

if __name__ == "__main__":

    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'mark2cure.settings')
    os.environ.setdefault('DJANGO_CONFIGURATION', 'Development')

    try:
        # from django.core.management import execute_from_command_line
        from configurations.management import execute_from_command_line

    except ImportError:
        try:
            import django
        except ImportError:
            raise ImportError(
                "Couldn't import Django. Are you sure it's installed and "
                "avilable on your PYTHONPATH environment variable? Did you "
                "forget to activate a viritual environment?"
            )
        raise
    execute_from_command_line(sys.argv)
