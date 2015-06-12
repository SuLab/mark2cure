import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "mark2cure.settings")
os.environ.setdefault("DJANGO_CONFIGURATION", "Production")

from configurations.wsgi import get_wsgi_application
application = get_wsgi_application()
