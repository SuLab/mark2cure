from allauth.account.utils import get_next_redirect_url
from allauth.account.views import SignupView

from ...task.models import Level
from django.utils import timezone


class SignupViewCustom(SignupView):

    def get_success_url(self):
        # After loggin them in, assign the first Level training so we know where to route them
        Level.objects.create(user=self.user, task_type=self.request.session.get('initial_training'), level=3, created=timezone.now())

        # Explicitly passed ?next= URL takes precedence
        ret = (
            get_next_redirect_url(
            self.request,
            self.redirect_field_name) or self.success_url)
        return ret

signup = SignupViewCustom.as_view()

