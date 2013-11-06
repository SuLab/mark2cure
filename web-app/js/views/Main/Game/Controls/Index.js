define(['marionette', 'templates', 'vent',
        //-- Models
        'models/User',
        'underscore.string'],
        function (Marionette, templates, vent,
                  User,
                  _s) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.main.game.controls,
    templateHelpers : function() { return this.options; },

    events : {
      'click button.done'       : 'submitAnnotations',
      'click button.navigate'   : 'nextDocument',
      'click .entity-sel button' : 'changeSelectionType'
    },

    initialize : function(options) {
      //-- This is the View w/ Control buttons and navigation
      //-- this.model == The Document
      //-- this.collection == The Collection of all Documents
      options.ann_list = _.uniq( this.model.get('annotations').pluck('text') ).sort();
      options.user = User;
      this.listenTo(this.model.get('annotations'), "add", this.reRender, this);
      this.listenTo(this.model.get('annotations'), "remove", this.reRender, this);
      this.listenTo(this.model, "change:matches", this.reRender, this);
      this.listenTo(this.options.user, "change:advance", this.reRender, this);
      this.listenTo(this.options.user, "change:sel_mode", this.reRender, this);
    },

    //
    //-- Events
    //
    reRender : function() {
      this.options.ann_list = _.uniq( this.model.get('annotations').pluck('text') ).sort();
      this.render();
    },

    submitAnnotations : function(evt) {
      evt.preventDefault();
      var c = true;

      if( this.options.ann_list.length ==  0) {
        c = confirm("Are you positive this document doesn't have any annotations?");
      }

      if(c) {
        this.model.save({'complete': true}, {wait: true});
        vent.trigger('network:refresh', {});
      }
    },

    nextDocument : function(evt) {
      evt.preventDefault();
      var index = this.collection.indexOf(this.model) + 1,
          index = (index < 0) ? this.collection.length-1 : index,
          index = (index == this.collection.length) ? 0 : index;

      this.options.model = this.model = this.collection.at(index);

      //-- Trigger the next view via actually changing url through router
      Backbone.history.navigate( '#/'+ this.model.id );
    },

    changeSelectionType : function(evt) {
      var type = $(evt.target).data('type');
      this.options.user.set('sel_mode', type);
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
