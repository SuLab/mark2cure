from allauth.account.utils import get_next_redirect_url
from allauth.account.views import SignupView

from ...task.models import Level
from django.utils import timezone


class SignupViewCustom(SignupView):

    def dispatch(self, request, *args, **kwargs):
        if not request.session.get('initial_training', False):
            return redirect('common:home')

        self.sociallogin = None
        data = request.session.get('socialaccount_sociallogin')
        if data:
            self.sociallogin = SocialLogin.deserialize(data)
        if not self.sociallogin:
            return HttpResponseRedirect(reverse('account_login'))
        return super(SignupView, self).dispatch(request, *args, **kwargs)

