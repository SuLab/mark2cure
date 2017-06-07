import requests
from requests import RequestException
from django.core.exceptions import PermissionDenied
from allauth.socialaccount.providers.base import ProviderException

from allauth.socialaccount import app_settings
from allauth.socialaccount.providers.base import AuthAction, AuthError
from allauth.socialaccount.models import SocialLogin, SocialToken
from allauth.socialaccount.helpers import (
    complete_social_login,
    render_authentication_error
)

from allauth.socialaccount.providers.oauth2.views import (
    OAuth2Adapter,
    OAuth2LoginView,
    OAuth2View
)
from allauth.socialaccount.providers.oauth2.client import (
    OAuth2Error,
)
from allauth.utils import build_absolute_uri, get_request_param


from .provider import ZooniverseProvider

class ZooniverseAdapter(OAuth2Adapter):
    provider_id = ZooniverseProvider.id

    # These are all production env endpoints
    access_token_url = 'https://www.zooniverse.org/oauth/token'
    authorize_url = 'https://panoptes.zooniverse.org/oauth/authorize'
    profile_url = 'https://www.zooniverse.org/api/me'
    redirect_uri_protocol = 'https'
    supports_state = False

    # User gets redirected to:
    # https://panoptes.zooniverse.org/oauth/authorize?response_type=token&client_id=17a0d&redirect_uri=https://localhost.getshopmatch.com:8000/accounts/zooniverse/login/callback
    #
    # Which then redirects to:
    # https://localhost.getshopmatch.com:8000/accounts/zooniverse/login/callback/?env=production#access_token=Md30PEUP9yisgCg6f0qI_aFaEkcww_iSlwemcM&token_type=bearer&expires_in=7200
    # (notice access_token + token_type=bearer + expires_in above - no refresh_token)

    # To refresh token: Open up the same login page (above) and parse the redirected URI (same domain) to get access_token

    # client.get_access_token(request.GET['code']) --->
    #   { 'access_token': '...', 'refresh_token': '...', 'expires_in': 1234 }

    settings = app_settings.PROVIDERS.get(provider_id, {})

    def complete_login(self, request, app, token, **kwargs):
        # Get extra info about the user (from the /api/me endpoint)
        headers = {
            "Authorization": "Bearer " + token.token,
            "Accept": "application/vnd.api+json; version=1" }
        if (self.headers): headers.update(self.headers)
        resp = requests.get(self.profile_url, headers=headers)
        resp_json = resp.json()

        return self.get_provider().sociallogin_from_response(
            request,
            resp_json['users'][0]
        )


class OAuth2CallbackView(OAuth2View):
    # Override the OAuth2CallbackView, since the Zooniverse authentication API doesn't support token_type=code
    def dispatch(self, request):
        if 'error' in request.GET or 'access_token' not in request.GET:
            # Distinguish cancel from error
            auth_error = request.GET.get('error', None)
            if auth_error == self.adapter.login_cancelled_error:
                error = AuthError.CANCELLED
            else:
                error = AuthError.UNKNOWN
            return render_authentication_error(
                request,
                self.adapter.provider_id,
                error=error)

        app = self.adapter.get_provider().get_app(self.request)
        client = self.get_client(request, app)
        try:
            # Parse the access token
            access_token = { 'access_token': request.GET['access_token'], 'expires_in': request.GET['expires_in'] }
            token = self.adapter.parse_token(access_token)
            token.app = app
            login = self.adapter.complete_login(request,
                                                app,
                                                token,
                                                response=access_token)
            login.token = token
            SocialLogin.stash_state(request)

            if self.adapter.supports_state:
                login.state = SocialLogin \
                    .verify_and_unstash_state(
                        request,
                        get_request_param(request, 'state'))
            else:
                login.state = SocialLogin.unstash_state(request)

            # Make sure the user is initiated with some training in the DB
            request.session['initial_training'] = 'zooniverse'

            return complete_social_login(request, login)
        except (PermissionDenied,
                OAuth2Error,
                RequestException,
                ProviderException) as e:
            return render_authentication_error(
                request,
                self.adapter.provider_id,
                exception=e)

oauth2_login = OAuth2LoginView.adapter_view(ZooniverseAdapter)
oauth2_callback = OAuth2CallbackView.adapter_view(ZooniverseAdapter)

