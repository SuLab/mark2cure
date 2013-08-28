define(['marionette', 
        'views/Main/Library/Item'], 
    function (Marionette, 
              LibItem) {
  'use strict';

  return Marionette.CollectionView.extend({
    itemView: LibItem,

  });
});
