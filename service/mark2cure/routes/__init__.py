from .. import factory

app = factory.create_app(__name__)

@app.route("/settings")
def settings():
    print "/ / / / / / /  setttingsss / / / / / /  / /"
