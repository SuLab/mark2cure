define(['marionette', 'templates', 'vent'],
    function (Marionette, templates, vent) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.main.game.results.index,
    tagName : "span",

    events : {
    },

    initialize : function(options) {
    },

    onRender : function() {
      var pop = this.model.get('parentDocument').get('popularity'),
          highlight = this.options.color_scale(pop[this.model.get('position')]);
      this.$el.css({'backgroundColor' : highlight});
    }
  });
});
