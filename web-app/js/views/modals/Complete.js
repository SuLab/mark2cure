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

    ui : {
      'email' : 'input'
    },

    events : {
      'click a.library'  : 'keepGoing',
      'blur input'       : 'submitEmail'
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
    },

    keepGoing : function(evt) {
      evt.preventDefault();

      vent.trigger('modal:close');
      Backbone.history.navigate( '#/library' );
      vent.trigger('navigate:library', {});
    }

  });
});
