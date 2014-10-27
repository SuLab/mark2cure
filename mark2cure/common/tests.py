

from django.contrib.auth.models import User
from django.test import TestCase, Client, LiveServerTestCase

from selenium.webdriver.firefox.webdriver import WebDriver


class SimpleTest(TestCase):

    def setUp(self):
        self.client = Client()

        user = User.objects.create_user('test-admin', 'test-admin@localhost.com', 'test-admin-password')
        user.is_staff = True
        user.save()


    def test_landing(self):
        '''
          Test that the landing page is online / rendering correctly
        '''
        response = self.client.get('/')
        self.assertTemplateUsed(response, 'landing/index.jade')


    def test_login(self):
        '''
          Make sure you need an account, redirect them if not authenticated
        '''
        response = self.client.get('/library/')
        self.assertRedirects(response, '/account/login/?next=/library/')


    def test_library_documents(self):
        '''
          Make sure documents are being served to the library view to iterate over
        '''
        logged_in = self.client.login(username = 'test-admin', password  ='test-admin-password')

        print logged_in

        response = self.client.get('/library')

        print response.status_code
        # check we've passed the polls to the template
        docs_in_context = response.content
        print "D O C S :"
        print docs_in_context
        # self.assertEquals(list(polls_in_context), [poll1, poll2])

        # check the poll names appear on the page
        # self.assertIn(poll1.question, response.content)
        # self.assertIn(poll2.question, response.content)


class MySeleniumTests(LiveServerTestCase):

    @classmethod
    def setUpClass(cls):
        cls.selenium = WebDriver()
        super(MySeleniumTests, cls).setUpClass()

    @classmethod
    def tearDownClass(cls):
        cls.selenium.quit()
        super(MySeleniumTests, cls).tearDownClass()

    def test_email_signup(self):
        '''
          Ensures that the user notify signup form and
          javascript is working and saving to the database
        '''

        demo_user_email = "demo@gmail.com"

        # Assert: No previous signup users
        signup_users = User.objects.filter(email = demo_user_email).all()
        self.assertEquals(len(signup_users), 0)

        # Submit the form
        self.selenium.get( self.live_server_url )
        email_input = self.selenium.find_element_by_name("email")
        email_input.send_keys(demo_user_email)
        self.selenium.find_element_by_name("email_notify").click()
        self.selenium.find_element_by_xpath('//button[@type="submit"]').click()

        # Assert: User added
        signup_users = User.objects.filter(email = demo_user_email).all()
        self.assertEquals(len(signup_users), 1)

        # Assert: email notify flag saved correctly
        signup_user = signup_users[0]
        self.assertEquals(signup_user.userprofile.email_notify, True)
