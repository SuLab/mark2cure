from django.utils import timezone
from .models import UserProfile

class ActiveUserMiddleware:

    def process_request(self, request):
        current_user = request.user
        if current_user.is_authenticated():
            UserProfile.objects.filter(user=current_user).update(last_seen=timezone.now())
