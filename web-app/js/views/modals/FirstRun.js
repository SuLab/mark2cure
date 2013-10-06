define(['marionette', 'templates', 'vent',
        //-- Models
        'models/User',],

        function ( Marionette, templates, vent,
                   User ) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.modals.first_run,
    className : 'modal-dialog',

    ui : {
      'experience'  : '#experience',
      'user_name'   : '#user-name',
    },

    events : {
      'click #myTab a'      : 'tabEvent',

      'change #experience'  : 'saveUserInfo',
      'blur #user-name'     : 'saveUserInfo',

      'click button'  : 'finish'
    },

    initialize : function(options) {
      this.model = options.user;
    },

    //
    //-- Events
    //
    tabEvent : function(evt) {
      evt.preventDefault()
      $(evt.target).tab('show')
    },

    finish : function(evt) {
      this.saveUserInfo();

      var doc = this.collection.at(0)
      Backbone.history.navigate( '#/document/'+ doc.id );
      vent.trigger('navigate:document', doc.id );
    },

    saveUserInfo : function(evt) {
      if(evt) { evt.preventDefault(); }
      this.model.set('experience',  Number(this.ui.experience.val()) );
      this.model.set('username',    this.ui.user_name.val()  );
      this.model.save();
    }

  });
});
