define(['backbone'],
    function(Backbone){
  'use strict';

  return Backbone.RelationalModel.extend({
    defaults: {
      text      : '',
      position  : 0,

      length  : null,
      start   : null,
      stop    : null,

      selected  : false,
      latest    : false,

      pos_x   : null,
      pos_y   : null,
      width   : null,
      height  : null,
    },

  });
});
