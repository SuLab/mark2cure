define(['vent'], 
    function (vent) {
  'use strict';

  return {
    showSettings : function() {
      vent.trigger('navigate:settings', {});
    },

    showFeedback : function() {
      vent.trigger('navigate:feedback', {});
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
