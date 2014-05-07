from django.contrib.auth.models import User


def get_mturk_account(worker_id):
    u, created = User.objects.get_or_create(username=worker_id)
    if created:
        u.set_password('')
        profile = u.profile
        profile.email_notify = False
        profile.mturk = True
        profile.save()
        u.save()

    return u


