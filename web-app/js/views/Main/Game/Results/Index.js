define(['marionette', 'templates', 'vent',
        //-- ETC
        'd3'], 

        function (Marionette, templates, vent,
                  Words) {
  'use strict';

  return Marionette.Layout.extend({
    template : templates.main.results.index,

    className : 'results-view',

    ui : {
    },

    events : {
    },

    initialize : function(options) {
      //-- View to show community consensus and self annotations
      //-- this.model == The Document Resource

      //-- i think these are in the this.model
      //-- this.collection == The Collection of user's annotations?
      console.log('results index init :: ', this);
    },

    onRender : function() {
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
    },

    //
    //-- Events
    //


    //
    //-- Utils
    //
    compare : function(user, truth) {
      var a = 0, b = 0, c = 0, d = 0, fpr, fnr;
      if( truth.length !== user.length ) return;
      for (var i = 0; i < truth.length; i++) {
        switch(user[i] + truth[i]) {
          case 0: ++d; break;
          case 2: ++a; break;
          case 1: user[i] ? ++c : ++b; break;
        }
      }
      if( a+b+c+d !== truth.length ) return;
      fpr = c/(a+c); fnr = b/(b+d);
      return {'true_pos': a, 'false_neg': b, 'false_pos': c, 'true_neg': d, 'false_pos_rate': fpr, 'false_neg_rate': fnr}
    }

  });
});
