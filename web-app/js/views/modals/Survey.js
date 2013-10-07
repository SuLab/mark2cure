define(['marionette', 'templates', 'vent',
        //-- Models
        'models/User',
        ],

        function ( Marionette, templates, vent,
                   //-- Models
                   User ) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.modals.survey,
    templateHelpers : function() { return this.options; },
    className : 'modal-dialog',

    events : {
      'click .close-btn'    : function(e) { e.preventDefault(); vent.trigger('modal:close'); },

      'change select' : 'changeInput'
    },

    initialize : function(options) {
      this.model = options.user;
      options['questions'] = [
        "I feel like I'm contributing to \"science\"",
        "I would do this regularly",
        "I felt comfortable doing these tasks",
        "I learned something doing this task?"]
    },

    //
    //-- Events
    //
    changeInput : function(evt) {
      evt.preventDefault();
      var val = Number($(evt.target).val()),
          q_index = $(evt.target).attr('name');

      this.model.set('feedback_'+q_index, val);
      this.model.save();
    }

  });
});
