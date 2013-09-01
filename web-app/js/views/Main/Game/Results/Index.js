define(['marionette', 'templates', 'vent',
        //-- ETC
        'd3'], 

        function (Marionette, templates, vent,
                  Words) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.main.game.results.index,
    tagName : 'p',
    className : 'paragraph results',

      // initialize : function(options) {
      //-- View to show community consensus and self annotations
      //-- this.model == The Document Resource

      //-- i think these are in the this.model
      //-- this.collection == The Collection of user's annotations?
      // console.log('results index init :: ', this);
    // },

    // onRender : function() {
      // var self = this,
      //     //------
      //     // HEATMAP CUTOFF %
      //     // -----
      //     community_threshold = .65;
      // this.options.max_pop = _.max( this.model.get('popularity') );
      // this.options.color_scale = d3.scale.linear()
      //                             .domain([0, this.options.max_pop])
      //                             .range(["white", "yellow"]);

      // if( this.model.get('complete') ) {
      //   var defined_correct = _.map(this.model.get('popularity'), function(i) { return i/self.options.max_pop >= community_threshold ? true : false; }),
      //       selected_index  = this.model.get('annotations').pluck('position'),
      //       user_correct    = _.map(this.model.get('popularity'), function(x, index) { return _.contains(selected_index, index); });
      //   this.options.score = this.compare(defined_correct, user_correct);
      //   this.options.ann_range = this.model.get('annotations').getRange();
      // }
      //
      //

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

      //   this.options.user.save({'first_run' : false});
      // }
    // },

    onClose : function() {
      //-- Incase they didn't close it before or are skipping
      $('.popover').hide();
    },


  });
});
