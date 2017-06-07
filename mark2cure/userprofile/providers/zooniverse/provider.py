from allauth.socialaccount.providers.base import ProviderAccount
from allauth.socialaccount.providers.oauth2.provider import OAuth2Provider
from allauth.socialaccount import providers
from allauth.account.models import EmailAddress


class ZooniverseAccount(ProviderAccount):
    def to_str(self):
        dflt = super(ZooniverseAccount, self).to_str()
        name = self.account.extra_data.get('display_name', dflt)
        return name


class ZooniverseProvider(OAuth2Provider):
    id = 'zooniverse'
    name = 'Zooniverse'
    account_class = ZooniverseAccount

    def extract_uid(self, data):
        # That's the user's unique identifier (his Zooniverse user ID)
        return data['id']

    def extract_email_addresses(self, data):
        # Extract the user's email address
        return [ EmailAddress(email=data['email'], verified=True, primary=True) ]

    def extract_common_fields(self, data):
        return dict(username=data.get('display_name'),
                    name=data.get('display_name'),
                    email=data.get('email'))


    def get_default_scope(self):
        scope = ['user', 'public']
        return scope


provider_classes = [ ZooniverseProvider ]
