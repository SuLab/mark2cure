from django.contrib.auth.forms import AuthenticationForm


def inject_signup_forms(request):
    return {} if request.user.is_authenticated() else {'login_form': AuthenticationForm()}
