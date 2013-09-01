define(['backbone', 'models/Word'], 
    
        function(Backbone, Word) {
  'use strict';

  return Backbone.Collection.extend({
    model   : Word,
    url     : '/api/v1/words',

    clear : function(attr) {
      return this.each(function(word) { word.set(attr, false); });
    },

  });
});
