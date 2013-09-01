define(['backbone', 'models/Word'], 
    
        function(Backbone, Word) {
  'use strict';

  return Backbone.Collection.extend({
    model   : Word,
    url     : '/api/v1/words',

    clearLatest : function() {
      return this.each(function(word) { word.set('latest', false); });
    },

    clearSelected : function() {
      return this.each(function(word) { word.set('selected', false); });
    }

  });
});
