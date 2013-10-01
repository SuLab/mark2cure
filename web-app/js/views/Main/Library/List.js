define(['marionette', 'templates', 'vent',
        'views/Main/Library/Item'], 
    function (Marionette, templates, vent,
              LibItem) {
  'use strict';

  return Marionette.CompositeView.extend({
    template : templates.main.library.list,
    itemView: LibItem,
    itemViewContainer : 'tbody'
  });
});
