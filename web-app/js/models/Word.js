define(['backbone'],
    function(Backbone){
  'use strict';

  return Backbone.RelationalModel.extend({
    defaults: {
      text      : '',

      length  : null,
      start   : null,
      stop    : null,

      latest    : false,
      selected  : false,
      neighbor    : false,
    },

  });
});
