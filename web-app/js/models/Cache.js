define(['backbone'],
    function(Backbone){
  'use strict';

  return Backbone.RelationalModel.extend({
    defaults: {
      keypress  : false,
    },

  });
});
