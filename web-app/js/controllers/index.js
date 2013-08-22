define(['vent'], 
    function (vent) {
  'use strict';

  return {
    showSettings : function() {
      vent.trigger('navigate:settings', {});
    },

    showInstructions : function() {
      vent.trigger('navigate:instructions', {});
    },

    showSurvey : function() {
      vent.trigger('navigate:survey', {});
    },

    sendMessage : function() {
      vent.trigger('navigate:message', {});
    },

    showLibrary : function() {
      vent.trigger('navigate:library', {});
    },

    setFilter : function(param) {
      vent.trigger('navigate', param || '');
    },
  };
});
