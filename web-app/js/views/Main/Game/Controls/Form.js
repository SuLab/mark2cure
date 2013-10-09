define(['marionette', 'templates', 'vent',
        //-- Models
        'models/User',
        'underscore.string'],
        function (Marionette, templates, vent,
                  User,
                  _s) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.main.game.form,
    templateHelpers : function() { return this.options; },

    ui : {
      'form' : 'form'
    },

    initialize : function(options) {
      //-- This is the View w/ Control buttons and navigation
      //-- this.model == The Document
      options.user = User;
    },

    onRender: function() {
      this.ui.form.submit();
    }

  });
});
