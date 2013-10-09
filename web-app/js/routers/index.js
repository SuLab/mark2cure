define(['marionette'],
    function(marionette) {
  'use strict';

  return marionette.AppRouter.extend({
    initialize: function() {
      this.bind('all', this._trackPageview);
    },

    _trackPageview: function() {
      var url;
      url = Backbone.history.getFragment();
      _gaq.push(['_trackPageview', "/" + url]);
    },

    appRoutes : {
      //-- Modal Pages
      'settings'      : 'showSettings',
      'instructions'  : 'showInstructions',

      'survey'    : 'showSurvey',
      'message'   : 'sendMessage',

      'library/:quest'  : 'showLibrary',
      'library/'         : 'showLibrary',
      'library'         : 'showLibrary',

      //-- Specific Document
      'document/:doc_id?assignmentId=:var1&hitId=:var2&workerId=:var3&turkSubmitTo=:var4' : 'showDocument',
      'document/:doc_id?assignmentId=:var1&hitId=:var2' : 'showDocument',
      'document/:doc_id' : 'showDocument',

      '*filter'         : 'showLibrary'
    }

  });
});
