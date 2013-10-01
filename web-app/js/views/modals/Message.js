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
    template : templates.modals.message,
    templateHelpers : function() { return this.options; },
    className : 'modal-dialog',

    ui : {
      'text'        : 'textarea#new_message',
      'char_count'  : 'span#char_count',
      'button'      : 'button'
    },

    events : {
      'click .close-btn' : function(e) { e.preventDefault(); vent.trigger('modal:close'); },

      'input textarea#new_message'  : 'calcLeft',
      'click button'                : 'sendMessage'
    },

    initialize : function(options) {
      this.model = options.user;
      options['char_limit'] = 320;
    },

    //
    //-- Events
    //
    calcLeft : function(evt) {
      evt.preventDefault();
      var remaining = this.options.char_limit - this.ui.text.val().length;
      this.ui.char_count.html( remaining );
      this.ui.button.attr("disabled", !(remaining >= 0) );
    },

    sendMessage : function(evt) {
      evt.preventDefault();
      var self = this;
      $.ajax({
        type: 'POST',
        url: '/api/v1/messages',
        dataType: 'json',
        data: {message : self.ui.text.val()},
        success: function () {
          self.ui.text.val('');
          self.close();
        }
      });
    }

  });
});
