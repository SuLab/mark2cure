define(['marionette', 'templates', 'vent',
        //-- Models & Collections
        'views/Main/Game/EntityTag/Words', 'views/Main/Game/Results/Index'],

        function (Marionette, templates, vent,
                  EntityTag, Results) {
  'use strict';

  return Marionette.Layout.extend({
    template : templates.main.game.index,
    templateHelpers : function() { return this.options; },

    regions : {
      game      : 'div.game',
    },

    ui : {
      'complete_alert'  : '.alert.alert-success',
      'navigate'        : 'button.navigate'
    },

    events : {
      'click button.done'   : 'submitAnnotations',
      'click button.navigate'    : 'nextDocument',
    },

    initialize : function(options) {
      console.log('Game Index :: ', this);
      //-- This is the Game view
      //-- this.model == The Document Resource
      //-- this.collection == The Collection of all Documents
      //-- options.user == The Currently Logged in User
    },

    onRender : function() {
      if( this.model.get('complete') ) {
        this.game.show( new Results({model: this.model}) );
      } else {
        this.game.show( new EntityTag(this.options) );
      }
    },

    //
    //-- Events
    //
    submitAnnotations : function(evt) {
      var self = this;
      evt.preventDefault();

      this.model.save({'complete': true});
      this.render();

      if( $('#network').length ) { vent.trigger('network:refresh', {}); }
      if( this.collection.completed().length == 5 ) {
        vent.trigger('navigate:show_complete');
      };

    },

    nextDocument : function(evt) {
      evt.preventDefault();
      var index = this.collection.indexOf(this.model) + 1,
          index = (index < 0) ? this.collection.length-1 : index,
          index = (index == this.collection.length) ? 0 : index;

      this.options.model = this.model = this.collection.at(index);
      Backbone.history.navigate( '#/'+ this.model.id );

      this.render();
    }

  });
});
