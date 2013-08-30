define(['marionette', 'templates', 'vent',
        'views/Main/Game/EntityTag/WordItem'], 
    function (Marionette, templates, vent,
              WordItem) {
  'use strict';

  return Marionette.CollectionView.extend({
    itemView: WordItem,
  });
});
