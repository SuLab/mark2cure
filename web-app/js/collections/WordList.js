define(['backbone', 'models/Word'], 
    
        function(Backbone, Word) {
  'use strict';

  return Backbone.Collection.extend({
    model   : Word,
    url     : '/api/v1/words',

    getSelected : function() {
      return this.where({selected : true});
    }

  });
});
