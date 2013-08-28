define(['marionette', 'templates', 'vent',
        'views/Main/Game/EntityTag/WordItem'], 
    function (Marionette, templates, vent,
              WordItem) {
  'use strict';

  return Marionette.CollectionView.extend({
    itemView: WordItem,

    initialize : function(options) {
      console.log('EntityTag Words :: ', this);
      this.collection = this.model.get('words');
    },

    // buildItemView: function(item, ItemView){
      // return new ItemView({ model : item,
                            // ann_range: this.options.ann_range });
    // }

  });
});
