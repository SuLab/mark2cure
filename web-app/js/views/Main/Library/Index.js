define(['marionette', 'templates', 'vent',
        //-- Views
        'views/Main/Library/List'],
    function (Marionette, templates, vent,
              //-- Views
              List) {
  'use strict';

  return Marionette.Layout.extend({
    template : templates.main.library.index,
    className : 'library-view',

    regions : {
      all   : '.library-list'
    },

    onRender : function() {
      this.all.show( new List({collection: this.collection}) );
    }

  });
});
