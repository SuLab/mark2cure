define(['backbone', 'models/Cache'], 
    
        function(Backbone, Cache) {
  'use strict';

  return Backbone.Collection.extend({
    model   : Cache,
    url     : '/api/v1/caches',

  });
});
