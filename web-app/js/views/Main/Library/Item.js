define(['marionette', 'templates', 'vent',
        'underscore.string'],

        function (Marionette, templates, vent,
                  _s) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.library.item,
    templateHelpers : function() { return this.options; },
    className: 'doc-item',

    events : {
      'click' : 'clickLibItem'
    },

    //
    //-- Events
    //
    clickLibItem : function(evt) {
      evt.preventDefault();
      Backbone.history.navigate( '#/'+ this.model.id );
      vent.trigger('navigate', this.model.id );
    }

  });
});
