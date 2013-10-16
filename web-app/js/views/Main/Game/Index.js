define(['marionette', 'templates', 'vent',
        //-- Views
        'views/Main/Game/EntityTag/Index', 'views/Main/Game/Results/Index', 'views/Main/Game/Controls/Index',
        'views/Main/Game/Controls/Form',
        //-- Models & Collections
        'models/User'],

        function (Marionette, templates, vent,
                  //-- Views
                  EntityTag, Results, Controls,
                  Form,
                  User) {
  'use strict';

  return Marionette.Layout.extend({
    template : templates.main.game.index,
    className : 'game-view',

    regions : {
      submission : 'div.submission',
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
      this.listenTo(this.model, "change:complete", this.reRender, this);
    },

    onRender : function() {
      this.reRender();
    },

    reRender : function() {
      this.controls.show( new Controls({model: this.model, collection: this.collection}) );

      if( this.model.get('complete') ) {
        //-- If you're showing the results and you're a turker, let amazon know to pay
        var highlighted = this.model.get('words').where({'selected': true});

        if( User.authenticated() &&
            User.get('mturk') &&
            User.get('assignment_id') &&
            highlighted.length
          ) {
          this.submission.show( new Form({model: this.model}) );
        }

        this.game.show( new Results({collection: this.model.get('words')}) );
      } else {
        if( User.get('assignment_id') == 'ASSIGNMENT_ID_NOT_AVAILABLE' ) {
          this.game.show( new Results({collection: this.model.get('words')}) );
        } else {
          this.game.show( new EntityTag({collection: this.model.get('words')}) );
        }
      }
    }

  });
});
