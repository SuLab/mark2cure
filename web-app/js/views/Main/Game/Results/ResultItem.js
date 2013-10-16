define(['marionette', 'templates', 'vent',
        'models/User'],
    function (Marionette, templates, vent,
             User) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.main.game.results.index,
    tagName : "span",

    onRender : function() {

      if(!User.get('mturk') || User.get('assignment_id')) {
        //-- Draw the community consensus
        var pop = this.model.get('parentDocument').get('popularity')[ this.model.collection.indexOf( this.model ) ];
        if(pop >= 1) {
          this.$el.html(this.model.get('text'));
          this.$el.addClass('neighbor');
          this.$el.css({'backgroundColor' : this.options.color_scale(pop)});
        }

        //-- Draw the user's annotations
        if( _.contains(this.options.ann_range, this.model.get('start')) ) {
          this.$el.addClass('selected');
          this.$el.addClass('neighbor');
        }
      }
    }

  });
});
