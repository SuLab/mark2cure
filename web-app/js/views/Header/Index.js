define(['marionette', 'templates', 'vent'],

        function (Marionette, templates, vent) {
  'use strict';

  return Marionette.Layout.extend({
    template : templates.header.index,
    templateHelpers : function() { return this.options; },
    className : 'navbar navbar-fixed-top',

    regions : {
      analytics : '#analytics-container'
    },

    events : {
      'click a.network'   : 'showNetwork',
      'click a.welcome'   : 'showWelcome',
    },

    initialize : function(options) {
      this.model = options.user;
      this.bindTo(this.collection,  'change:complete', this.render, this);
      this.bindTo(this.model,       'all', this.render, this);
    },

    //
    //-- Events
    //
    showNetwork : function(evt) {
      evt.preventDefault();
      vent.trigger('navigate:analytics', {toggle: true});
    },

    showWelcome : function(evt) {
      evt.preventDefault();
      vent.trigger('modal:welcome', {});
    }

  });
});
