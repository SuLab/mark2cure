import raven
import os

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
PROJECT_PATH = os.path.realpath(os.path.join(os.path.dirname(__file__), os.path.pardir))

ADMINS = ()
MANAGERS = ADMINS
SITE_ID = 1
INTERNAL_IPS = ('127.0.0.1',)

# Application definition
INSTALLED_APPS = (
    'sslserver',  # HTTPS local development server
    'grappelli',
    'django.contrib.admin',

    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sites',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',

    'django_comments',

    'django.contrib.flatpages',
    'django.contrib.sitemaps',
    'raven.contrib.django.raven_compat',

    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    'allauth.socialaccount.providers.google',
    'mark2cure.userprofile.providers.zooniverse',

    'djrill',
    # 'djcelery',
    'robots',
    'corsheaders',

    'django_extensions',
    # 'django_nose',
    'rest_framework',
    # 'debug_toolbar',

    # Mark2Cure specific apps
    'mark2cure.document',

    'mark2cure.userprofile',
    'mark2cure.team',
    'mark2cure.instructions',
    'mark2cure.training',

    'mark2cure.task',
    'mark2cure.task.entity_recognition',
    'mark2cure.task.relation',

    'mark2cure.score',
    'mark2cure.common',
    'mark2cure.talk',
    'mark2cure.api',
    'mark2cure.download',
    'mark2cure.control',
    'mark2cure.analysis',

    'django.contrib.humanize',
    'widget_tweaks',
    'storages',
    'gunicorn'
)

LOGGING = {
    'version': 1,
    'disable_existing_loggers': True,
    'root': {
        'level': 'WARNING',
        'handlers': ['sentry'],
    },
    'formatters': {
        'verbose': {
            'format': '%(levelname)s %(asctime)s %(module)s '
                      '%(process)d %(thread)d %(message)s'
        },
        # 'django.server': {
        #     '()': 'django.utils.log.ServerFormatter',
        #     'format': '[%(server_time)s] %(message)s',
        # }
    },
    'handlers': {
        'sentry': {
            'level': 'ERROR',
            'class': 'raven.contrib.django.raven_compat.handlers.SentryHandler',
            'tags': {'custom-tag': 'x'},
        },
        'console': {
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'verbose'
        },
        # 'django.server': {
        #     'level': 'INFO',
        #     'class': 'logging.StreamHandler',
        #     'formatter': 'django.server',
        # },
    },
    'loggers': {
        'django.db.backends': {
            'level': 'ERROR',
            'handlers': ['console'],
            'propagate': False,
        },
        'raven': {
            'level': 'DEBUG',
            'handlers': ['console'],
            'propagate': False,
        },
        'sentry.errors': {
            'level': 'DEBUG',
            'handlers': ['console'],
            'propagate': False,
        },
        # 'django.server': {
        #     'handlers': ['django.server', 'console'],
        #     'level': 'INFO',
        #     'propagate': False,
        # }
    },
}

AUTHENTICATION_BACKENDS = (
    'django.contrib.auth.backends.ModelBackend',

    # `allauth` specific authentication methods, such as login by e-mail
    'allauth.account.auth_backends.AuthenticationBackend',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.cache.UpdateCacheMiddleware',

    'django.contrib.sessions.middleware.SessionMiddleware',
    'corsheaders.middleware.CorsMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    # 'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'mark2cure.userprofile.activeuser_middleware.ActiveUserMiddleware',

    'raven.contrib.django.raven_compat.middleware.Sentry404CatchMiddleware',

    'django.middleware.cache.FetchFromCacheMiddleware',

    # 'account.middleware.LocaleMiddleware',
    # 'account.middleware.TimezoneMiddleware',
)

ROOT_URLCONF = 'mark2cure.urls'
WSGI_APPLICATION = 'mark2cure.wsgi.application'

PASSWORD_HASHERS = (
    'django.contrib.auth.hashers.BCryptSHA256PasswordHasher',
    'django.contrib.auth.hashers.BCryptPasswordHasher',
    'django.contrib.auth.hashers.PBKDF2PasswordHasher',
    'django.contrib.auth.hashers.PBKDF2SHA1PasswordHasher',
)

LANGUAGE_CODE = 'en-us'
TIME_ZONE = 'America/Los_Angeles'
USE_I18N = True
USE_L10N = True
USE_TZ = True

# Auto signup a user after social login (using Zooniverse) - so the user won't have to go through our sign-up form after initial login
AUTO_SIGNUP = True

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            PROJECT_PATH + '/templates/',
        ],
        'OPTIONS': {
            'debug': True,
            'context_processors': [
                'django.template.context_processors.request',

                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.i18n',
                'django.template.context_processors.media',
                'django.template.context_processors.static',
                'django.template.context_processors.tz',
                'django.contrib.messages.context_processors.messages',
                # App Specific
                'mark2cure.common.context_processors.support_form',
            ],
            'loaders': [
                ('pyjade.ext.django.Loader', (
                    'django.template.loaders.filesystem.Loader',
                    'django.template.loaders.app_directories.Loader',
                )),
            ],
            'builtins': ['pyjade.ext.django.templatetags'],
        },
    },
]


MEDIA_URL = 'media/'
MEDIA_ROOT = 'media/'

STATIC_URL = '/static/'
STATIC_ROOT = '/static/'

STATICFILES_DIRS = (
    PROJECT_PATH + '/static',
    'static'
)

FIXTURE_DIRS = (
    PROJECT_PATH + '/fixtures/',
)

DEFAULT_FILE_STORAGE = 'storages.backends.s3boto.S3BotoStorage'
AWS_STORAGE_BUCKET_NAME = 'mark2cure'
AWS_QUERYSTRING_AUTH = False

STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
    'django.contrib.staticfiles.finders.FileSystemFinder',
)

GRAPPELLI_AUTOCOMPLETE_SEARCH_FIELDS = True

# Admin/Control settings
ADMIN_MEDIA_PREFIX = STATIC_URL + 'admin/'
SOUTH_TESTS_MIGRATE = False
REUSE_DB = 1

EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'

# User settings
AUTH_PROFILE_MODULE = 'userprofile.UserProfile'
LOGIN_URL = '/accounts/login/'
LOGOUT_URL = '/accounts/logout/'
LOGIN_REDIRECT_URL = '/'

REST_FRAMEWORK = {
    'DEFAULT_AUTHENTICATION_CLASSES': (
        'rest_framework.authentication.BasicAuthentication',
        'rest_framework.authentication.SessionAuthentication',
    )
}

# User is online if they've been last seen 5min ago
USER_ONLINE_TIMEOUT = 300

ROBOTS_USE_SITEMAP = True

# User account management
ACCOUNT_AUTHENTICATION_METHOD = 'username_email'
ACCOUNT_EMAIL_REQUIRED = True
ACCOUNT_EMAIL_CONFIRMATION_ANONYMOUS_REDIRECT_URL = ' (TODO) '
ACCOUNT_EMAIL_CONFIRMATION_EXPIRE_DAYS = 30
ACCOUNT_EMAIL_REQUIRED = True
ACCOUNT_SIGNUP_PASSWORD_VERIFICATION = True
ACCOUNT_UNIQUE_EMAIL = True
ACCOUNT_DEFAULT_HTTP_PROTOCOL = 'https'
ACCOUNT_USERNAME_MIN_LENGTH = 5
ACCOUNT_USERNAME_BLACKLIST = ['admin', 'owner', 'user']
ACCOUNT_USERNAME_REQUIRED = True
ACCOUNT_PASSWORD_MIN_LENGTH = 6
ACCOUNT_LOGIN_ON_EMAIL_CONFIRMATION = True
ACCOUNT_LOGIN_ON_PASSWORD_RESET = True
ACCOUNT_SESSION_REMEMBER = True
ACCOUNT_TEMPLATE_EXTENSION = 'jade'
SOCIALACCOUNT_QUERY_EMAIL = ACCOUNT_EMAIL_REQUIRED
SOCIALACCOUNT_PROVIDERS = {'google': {
    'SCOPE': ['profile', 'email'],
    'AUTH_PARAMS': {'access_type': 'online'}
}}

RAVEN_CONFIG = {'dsn': '', 'release': raven.fetch_git_sha(BASE_DIR)}
ENTREZ_TERMS = ['chordoma', 'breast cancer', 'diabetes']
ENTREZ_MAX_COUNT = 10

NCBO_API_KEY = '5dbf6de8-fda5-4ef7-a3a8-bcea293f3715'
STOP_WORDS = 'protein,gene,disease,disorder,syndrome,chromosome,receptor,cell,\
        orphan,thumb,center,with,involved,image,type,known,encoded,this,both,\
        human,second,near,observed,from,family,width,name,caption,state,\
        structure,MEROPS,Pfam,domain,Symbol,SMART,crystal,analogue,\
        protein family,SCOP,InterPro,EC number,Name,group,related,then,Some,\
        form,http,abstract,content,liter,levels,enzyme,drugs,into,slow,\
        intermediate,bound,Citation,when,down,After'

POSTIVE_FLATTER = ['Bravo', 'Wow', 'Super', 'Terrific', 'Cool', 'Amazing', 'Superb', 'Brilliant', 'Fantastic', 'Fabulous', 'You\'re a Champion', 'Well done', 'You rock', 'Great job', 'Tip top', 'Good thinking', 'Keep it up', 'Way to go', 'Right on', 'Top stuff', 'Take a bow', 'Unreal', 'Impressed', 'Great stuff', 'Awesome', 'Nice going', 'Very Creative', 'Thank You', 'Beautiful', 'Very proud', 'Good for you', 'Give me five', 'You make me happy', 'A+', 'A++', 'AA+', 'Fab', 'Rad', 'A+++', 'AAA+', 'A-OK', 'Best', 'Cool', 'Deal', 'Fast', 'Fine', 'Item', 'Nice', 'Safe', 'Thx!', 'WOW!', 'Prime', 'Solid', 'Super', 'Sweet', 'Thanx', 'Whoa!', 'Groovy', 'Honest', 'Speedy', 'Superb', 'Sweeet', 'Thanks', 'Zowie!', 'Amazing', 'Awesome', 'Quality', 'Service', 'Sweeeet', 'Glorious', 'Stunning', 'Superior', 'The Best', 'The Bomb', 'Thrilled', 'Way Cool', 'Brilliant', 'Competent', 'Delighted', 'Excellent', 'Exquisite', 'Marvelous', 'Overjoyed', 'Satisfied', 'Thank You', 'Top Notch', 'Unrivaled', 'Wonderful', 'A Home Run', 'Astounding', 'Delightful', 'Impressive', 'Incredible', 'Super Cool', 'Super Fast', 'Supersonic', 'Astonishing', 'Fascinating', 'Interesting', 'Magnificent', 'No Problems', 'Outstanding', 'Splendorous', 'Trustworthy', 'Unsurpassed', 'Wicked Cool', 'Breathtaking', 'Looking Good', 'Overwhelming!', 'Unbelievable!', 'Awe Inspiring', 'Lickety Split', 'Splendiferous', 'Thanks A Ton!', 'Extremely Cool', 'Satisfied 100%', 'Extremely Happy', 'Great Condition', 'Above And Beyond', 'State Of The Art', 'Thanks A Million!', 'Unbelievably Cool', 'Expertly Described', 'Extremely Satisfied', 'Great Communication', 'Greatly Appreciated', 'Beyond My Wildest Dreams', 'Supercalifragilisticexpialidocious', 'Thank You! Thank You! Thank You!']
SUPPORT_FLATTER = ['You can do it', 'Nice Try', 'Don\'t give up', 'Every bit counts', 'Thank you', 'Keep going', 'You can do better than that']
ENTITY_RECOGNITION_K = 15
RELATION_K = 15
ENTITY_RECOGNITION_DOC_POINTS = 1000

RELATION_REL_POINTS = 75
RELATION_DOC_POINTS = 1000

# Email settings management
DEFAULT_FROM_EMAIL = 'Mark2Cure <contact@mark2cure.org>'
SERVER_EMAIL = DEFAULT_FROM_EMAIL
EMAIL_BACKEND = "djrill.mail.backends.djrill.DjrillBackend"

# CELERY_DEFAULT_QUEUE = 'default'
# CELERY_DEFAULT_EXCHANGE_TYPE = 'direct'
# CELERY_DEFAULT_ROUTING_KEY = 'default'
# CELERY_DEFAULT_DELIVERY_MODE = 'persistent'
# # (TODO) Look into setting this up, same cert as HTTPS?
# # BROKER_USE_SSL // http://celery.readthedocs.io/en/latest/userguide/security.html#message-signing
# CELERY_RESULT_BACKEND = 'djcelery.backends.database:DatabaseBackend'
# CELERY_TIMEZONE = 'America/Los_Angeles'
# CELERYD_CONCURRENCY = 2
# CELERYD_PREFETCH_MULTIPLIER = 2  # limit for the number of tasks (messages) a worker can reserve
# CELERY_TASK_RESULT_EXPIRES = 3600  # If results are sent back, that they expire after one hour if not picked up/acknowledged.
# CELERY_MAX_CACHED_RESULTS = 100
# CELERY_TRACK_STARTED = True
# CELERY_TASK_SERIALIZER = 'json'
# CELERY_EVENT_SERIALIZER = 'json'
# CELERY_ACCEPT_CONTENT = ['application/json']
# # These are just global standards
# CELERYD_TASK_TIME_LIMIT = 60 * 30
# CELERYD_TASK_SOFT_TIME_LIMIT = 60 * 5
# CELERY_STORE_ERRORS_EVEN_IF_IGNORED = True
# CELERY_SEND_EVENTS = True  # Send events so the worker can be monitored by tools like celerymon.
# CELERY_SEND_TASK_SENT_EVENT = True
# CELERY_EVENT_QUEUE_TTL = 10 * 60  # Message expiry time in seconds (int/float) for when messages sent to a monitor clients event queue is deleted
# CELERY_RESULT_PERSISTENT = True
# CELERY_ALWAYS_EAGER = False


try:
    from .local_settings import *  # noqa
except ImportError as e:
    if "local_settings" not in str(e):
        raise e

