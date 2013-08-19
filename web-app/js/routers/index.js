define(['marionette'], 
    function(marionette) {
  'use strict';

  return marionette.AppRouter.extend({

    appRoutes:{
      //-- Modal Pages
      'settings' : 'showSettings',
      'feedback' : 'showFeedback',
      'message'  : 'sendMessage',

      //-- Pages
      'library'  : 'showLibrary',

      //-- Etc
      '*filter'  : 'setFilter'
    }

  });
});
