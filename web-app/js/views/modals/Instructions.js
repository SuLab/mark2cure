define(['marionette', 'templates', 'vent'], 

        function ( Marionette, templates, vent ) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.modals.instructions,
    className : 'modal-dialog',

    events : {
      'click .close-btn' : function(e) { e.preventDefault(); vent.trigger('modal:close'); },
    }

  });
});
