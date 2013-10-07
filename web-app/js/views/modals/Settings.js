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
      'click button' : 'saveUser'
    },

    initialize : function(options) {
      this.model = options.user;
      this.listenTo(this.model, 'sync', this.render);
    },

    //
    //-- Events
    //
    saveUser : function(evt) {
      evt.preventDefault();
      this.model.set({'experience' : Number(this.ui.experience.val()) });
      this.model.set({'username' : this.ui.user_name.val()});
      this.model.set({'email' : this.ui.email.val()});
      this.model.set({'advance' : !this.model.get('advance')});
      this.model.save(null, {success: function() {
        vent.trigger('modal:close');
      }});
    },


  });
});
