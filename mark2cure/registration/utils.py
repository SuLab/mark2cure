try:
    from urllib.parse import urlunparse
except ImportError:
    from urlparse import urlunparse


def get_local_host(request):
    scheme = 'http' + ('s' if request.is_secure() else '')
    return url(scheme=scheme, host=request.get_host())


def url(scheme='', host='', path='', params='', query='', fragment=''):
    return urlunparse((scheme, host, path, params, query, fragment))

