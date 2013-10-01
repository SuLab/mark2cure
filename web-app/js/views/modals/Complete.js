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
    template : templates.modals.complete,
    className : 'modal-dialog',

    ui : {
      'email' : 'input',
      'email_btn' : '.btn'
    },

    events : {
      'click .close-btn' : function(e) { e.preventDefault(); vent.trigger('modal:close'); },
      'click .close-link' : function(e) { e.preventDefault(); vent.trigger('modal:close'); },

      'blur input'       : 'submitEmail',
      'click .btn.btn-success' : function(e) { e.preventDefault(); this.ui.email_btn.removeClass('btn-success').addClass('btn-info'); }
    },

    initialize : function(options) {
      this.model = options.user;
    },

    //
    //-- Events
    //
    submitEmail : function(evt) {
      evt.preventDefault();
      this.model.save({'email' : this.ui.email.val().trim()});
    }

  });
});
