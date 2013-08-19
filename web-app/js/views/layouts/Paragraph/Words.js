define(['marionette', 'views/layouts/Paragraph/WordItem'], 
    function (Marionette, WordItem) {
  'use strict';

  return Marionette.CollectionView.extend({
    itemView: WordItem,

    initialize : function(options) {
      this.collection = this.model.get('words');
    },

    buildItemView: function(item, ItemView){
      return new ItemView({ model : item,
                            max_pop : this.options.max_pop,
                            color_scale: this.options.color_scale,
                            ann_range: this.options.ann_range });
    }

  });
});
