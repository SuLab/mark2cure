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

    appRoutes:{
      //-- Modal Pages
      'settings'      : 'showSettings',
      'instructions'  : 'showInstructions',

      'survey'    : 'showSurvey',
      'message'   : 'sendMessage',

      //-- Pages
      'library'         : 'showLibrary',
      'library/:quest'  : 'showLibrary',

      //-- Etc
      '*filter'  : 'setFilter'
    }

  });
});
