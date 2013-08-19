define(['marionette', 'templates', 'vent',
        //-- Models
        'models/User',
        //-- libs
        ], 

        function ( Marionette, templates, vent,
                   //-- Models
                   User ) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.modals.settings,
    templateHelpers : function() { return this.options; },

    ui : {
      'experience'  : '#user-experience',
      'user_name'   : '#user-username',
      'email'       : '#user-email',
    },

    events : {
      'click .close-btn' : function(e) { e.preventDefault(); vent.trigger('modal:close'); },

      'change #user-experience'   : 'saveUserInfo',
      'blur #user-username'       : 'saveUserInfo',
      'blur #user-email'          : 'saveUserInfo'
    },

    initialize : function(options) {
      this.model = options.user;
    },

    //
    //-- Events
    //
    saveUserInfo : function(evt) {
      evt.preventDefault();
      this.model.set('experience',  Number(this.ui.experience.val()) );
      this.model.set('username',    this.ui.user_name.val()  );
      this.model.set('email',       this.ui.email.val()  );
      this.model.save();
    }

  });
});
