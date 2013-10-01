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
    className : 'modal-dialog',

    ui : {
      'experience'  : '#user-experience',
      'user_name'   : '#user-username',
      'email'       : '#user-email',
    },

    events : {
      'click .close-btn' : function(e) { e.preventDefault(); vent.trigger('modal:close'); },

      'change #user-experience'   : 'saveExperience',
      'blur #user-username'       : 'saveUserName',
      'blur #user-email'          : 'saveEmail',
      'change #advance'           : 'changeAdvanceSettings'
    },

    initialize : function(options) {
      this.model = options.user;
      this.listenTo(this.model, 'sync', this.render);
    },

    //
    //-- Events
    //
    saveExperience : function(evt) {
      evt.preventDefault();
      this.model.save({'experience' : Number(this.ui.experience.val()) });
    },

    saveUserName : function(evt) {
      evt.preventDefault();
      this.model.save({'username' : this.ui.user_name.val()});
    },

    saveEmail : function(evt) {
      evt.preventDefault();
      this.model.save({'email' : this.ui.email.val()});
    },

    changeAdvanceSettings : function() {
      this.model.save({'advance' : !this.model.get('advance')});
    }

  });
});
