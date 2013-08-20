define(['marionette', 'templates', 'vent',
        //-- Models
        'models/User',
        //-- libs
        'wizard'], 

        function ( Marionette, templates, vent,
                   User, wizard ) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.modals.first_run,

    ui : {
      'wizard'      : '#MyWizard',
      'experience'  : '#experience',
      'user_name'   : '#user-name',
    },

    events : {
      'change #experience'  : 'saveUserInfo',
      'blur #user-name'     : 'saveUserInfo',
      'keydown #user-name'  : 'finish'
    },

    initialize : function(options) {
      this.model = options.user;
    },

    onRender : function() {
      var self = this;
      this.ui.wizard
        .on('finished', function (e) {

          var doc = self.collection.at(0)
          Backbone.history.navigate( '#/'+ doc.id );
          vent.trigger('navigate', doc.id );

          //-- Completed training, hide background
          self.model.save({'first_run' : false});
          self.close();
        });
    },

    //
    //-- Events
    //
    finish : function(evt) {
      if(evt.keyCode === 13) {
        this.saveUserInfo();

        var doc = this.collection.at(0)
        Backbone.history.navigate( '#/'+ doc.id );
        vent.trigger('navigate', doc.id );

        //-- Completed training, hide background
        this.model.save({'first_run' : false});
        this.close();
      }
    },

    saveUserInfo : function(evt) {
      if(evt) { evt.preventDefault(); }
      this.model.set('experience',  Number(this.ui.experience.val()) );
      this.model.set('username',    this.ui.user_name.val()  );
      this.model.save();
    }

  });
});
