define(['marionette', 'bootstrap'], 
        function(marionette, bootstrap) {
  'use strict'

  return marionette.Region.extend({
    open: function(view) {
      this.$el.hide();
      this.$el.html(view.el);
      this.$el.slideDown('slow');
    }

  });
});
