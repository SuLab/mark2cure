define(['marionette', 'templates', 'vent',
        //-- Models
        'models/User',
        //-- ETC
        'moment'],
    function (Marionette, templates, vent,
              User) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.main.relationgame.objectitem,
    tagName : 'li',

  });
});
