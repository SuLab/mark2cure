define(['marionette', 'templates', 'vent',
        'views/Main/Game/Results/ResultItem',
        //-- ETC
        'd3'], 

        function (Marionette, templates, vent,
                  ResultItem) {
  'use strict';

  return Marionette.CollectionView.extend({
    itemView: ResultItem,
    tagName : 'p',
    className : 'paragraph results',

    initialize : function(options) {
      //-- View to show community consensus and self annotations
      //-- this.collection == The List of words Resource

      //------
      // HEATMAP CUTOFF %
      // -----
      var community_threshold = .65;
      this.options.popularity = _.max( this.collection.parentDocument.get('popularity') );
      this.options.color_scale = d3.scale.linear()
                                  .domain([0, this.options.popularity])
                                  .range(["white", "yellow"]);

      //   var defined_correct = _.map(this.model.get('popularity'), function(i) { return i/self.options.max_pop >= community_threshold ? true : false; }),
      //       selected_index  = this.model.get('annotations').pluck('position'),
      //       user_correct    = _.map(this.model.get('popularity'), function(x, index) { return _.contains(selected_index, index); });
      //   this.options.score = this.compare(defined_correct, user_correct);
      //   this.options.ann_range = this.model.get('annotations').getRange();

      // if( this.collection.completed().length == 1 ) {
      //   this.ui.navigate.popover({  title   : 'Congrats! You annotated your first document!',
      //                               content : templates.snippets.paragraph_info({}),
      //                               html    : true,
      //                               trigger : 'manual',
      //                               placement : 'left',
      //                               container : 'body' });
      //   this.ui.navigate.popover('show');

      //   $('button.explain-network').click(function(e) {
      //     self.ui.navigate.popover('hide');
      //     vent.trigger('navigate:analytics', {toggle: true, explain: true});
      //   });

      // }
    },

    onClose : function() {
      //-- Incase they didn't close it before or are skipping
      $('.popover').hide();
    },

    buildItemView: function(item, ItemViewType) {
      // build the final list of options for the item view type
      var options = _.extend({model: item}, { color_scale : this.options.color_scale });
      // create the item view instance
      var view = new ItemViewType(options);
      return view;
    }

  });
});
