define(['marionette', 
    'views/layouts/Library/Item'], 
    function (Marionette, 
              LibItem) {
  'use strict';

  return Marionette.CollectionView.extend({
    itemView: LibItem,

  });
});
