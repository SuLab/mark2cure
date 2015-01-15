from .forms import SupportMessageForm

def support_form(request):
    return {} if request.user.is_authenticated() else {'support_form': SupportMessageForm()}
