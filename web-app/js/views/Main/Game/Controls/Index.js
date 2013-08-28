define(['marionette', 'templates', 'vent'],
        function (Marionette, templates, vent) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.main.game.controls,

    events : {
      'click button.done'       : 'submitAnnotations',
      'click button.navigate'   : 'nextDocument'
    },

    initialize : function() {
      //-- This is the View w/ Control buttons and navigation
      //-- this.model == The Document
      //-- this.collection == The Collection of all Documents
    },

    //
    //-- Events
    //
    submitAnnotations : function(evt) {
      evt.preventDefault();
      this.model.save({'complete': true});

      // if( $('#network').length ) { vent.trigger('network:refresh', {}); }
      // if( this.collection.completed().length == 5 ) {
      //   vent.trigger('navigate:show_complete');
      // };
    },

    nextDocument : function(evt) {
      evt.preventDefault();
      var index = this.collection.indexOf(this.model) + 1,
          index = (index < 0) ? this.collection.length-1 : index,
          index = (index == this.collection.length) ? 0 : index;

      this.options.model = this.model = this.collection.at(index);

      //-- Trigger the next view via actually changing url through router
      Backbone.history.navigate( '#/'+ this.model.id );
    }

    //
    //-- Utils
    //
    // compare : function(user, truth) {
    //   var a = 0, b = 0, c = 0, d = 0, fpr, fnr;
    //   if( truth.length !== user.length ) return;
    //   for (var i = 0; i < truth.length; i++) {
    //     switch(user[i] + truth[i]) {
    //       case 0: ++d; break;
    //       case 2: ++a; break;
    //       case 1: user[i] ? ++c : ++b; break;
    //     }
    //   }
    //   if( a+b+c+d !== truth.length ) return;
    //   fpr = c/(a+c); fnr = b/(b+d);
    //   return {'true_pos': a, 'false_neg': b, 'false_pos': c, 'true_neg': d, 'false_pos_rate': fpr, 'false_neg_rate': fnr}
    // }

  });
});
