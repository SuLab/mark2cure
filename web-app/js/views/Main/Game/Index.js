define(['marionette', 'templates', 'vent',
        //-- Views
        'views/Main/Game/EntityTag/Words',
        'views/Main/Game/Results/Index',
        'views/Main/Game/Controls/Index'],

        function (Marionette, templates, vent,
                  //-- Views
                  EntityTag, 
                  Results, 
                  Controls) {
  'use strict';

  return Marionette.Layout.extend({
    template : templates.main.game.index,
    templateHelpers : function() { return this.options; },

    regions : {
      game      : 'div.game',
      controls  : 'div.controls'
    },

    ui : {
      'complete_alert'  : '.alert.alert-success',
      'navigate'        : 'button.navigate'
    },

    initialize : function(options) {
      //-- This is the Game view
      //-- this.model == The Document
      //-- this.collection == The Collection of all Documents
      //-- options.user == The Currently Logged in User

      if(!this.model.get('words').length) { this.model.parseText(); }
      this.listenTo(this.model, "change:complete", this.render, this);
    },

    onRender : function() {
      console.log('game index rerendzzz', this);
      this.controls.show( new Controls({model: this.model, collection: this.collection}) );

      if( this.model.get('complete') ) {
        this.game.show( new Results({model: this.model}) );
      } else {
        this.game.show( new EntityTag({model: this.model}) );
      }
    }

    //
    //-- Events
    //


  });
});
