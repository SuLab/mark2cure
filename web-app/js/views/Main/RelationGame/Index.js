define(['marionette', 'templates', 'vent',
        'views/Main/RelationGame/ObjectList'],
        function (Marionette, templates, vent,
                  ObjectList) {
  'use strict';

  return Marionette.Layout.extend({
    template : templates.main.relationgame.index,
    className : 'relationship-game-view row',

    regions : {
      'objects' : '.objects'
    },

    ui : {
    },

    initialize : function(options) {
      //-- This is the Game view
      //-- this.model == The Document
      //-- this.collection == The Collection of all Documents
      //-- options.user == The Currently Logged in User
      if(!this.model.get('words').length) { this.model.parseText(); }
    },

    onRender : function() {
      this.objects.show( new ObjectList({collection: this.model.get('annotations')}) );
    }


  });
});
