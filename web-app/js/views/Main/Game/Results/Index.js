define(['marionette', 'templates', 'vent',
        'views/Main/Game/Results/ResultItem'],

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
      var community_threshold = .65;
      // -----
      this.options.popularity = _.max( this.collection.parentDocument.get('popularity') );
      this.options.color_scale = d3.scale.linear()
                                  .domain([0, this.options.popularity])
                                  .range(["white", "yellow"]);
      this.options.ann_range = this.collection.parentDocument.get('annotations').getRange();
    },

    onRender : function() {
      if( this.collection.parentDocument.collection.completed().length === 1 ) {
        var $navigate = $('button.navigate');
        $navigate.popover({  title   : 'Congrats! You annotated your first document!',
                                    content : templates.snippets.paragraph_info({}),
                                    html    : true,
                                    trigger : 'manual',
                                    placement : 'left',
                                    container : 'body' });
       $navigate.popover('show');
      $('.popover').css({'position': 'absolute', top: 160, left: ($('body').width()/2)+80 });

        $('button.explain-network').click(function(e) {
          $('.popover').hide();
          vent.trigger('navigate:analytics', {toggle: true, explain: true});
        });
      }
    },

    onClose : function() {
      //-- Incase they didn't close it before or are skipping
      $('.popover').hide();
    },

    buildItemView: function(item, ItemViewType) {
      // build the final list of options for the item view type
      var options = _.extend({model: item}, { color_scale : this.options.color_scale,
                                              ann_range : this.options.ann_range });
      // create the item view instance
      var view = new ItemViewType(options);
      return view;
    }

  });
});
