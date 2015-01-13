try:
    from urllib.parse import parse_qs, urlencode, urljoin, urlunparse
except ImportError:
    from urllib import urlencode
    from urlparse import parse_qs, urljoin, urlunparse

from django.conf import settings


def get_local_host(request):
    scheme = 'http' + ('s' if request.is_secure() else '')
    return url(scheme=scheme, host=request.get_host())


def url(scheme='', host='', path='', params='', query='', fragment=''):
    return urlunparse((scheme, host, path, params, query, fragment))

