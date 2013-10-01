define(['marionette', 'templates', 'vent',
        'underscore.string'],

        function (Marionette, templates, vent,
                  _s) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.main.library.item,
    tagName: 'tr',

    events : {
      // 'click' : 'clickLibItem'
    },

    initialize : function() {
      this.listenTo(this.model, "change", this.render, this);
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
